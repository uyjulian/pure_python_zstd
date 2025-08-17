# Some numerical data is initialized as -1 even when it doesn't need initialization to help the JIT infer types

from typing import Union

# Huffman decoding table
class HDT:
	# initial bits
	b: int = 0
	# symbols
	s: memoryview = memoryview(bytearray())
	# num bits
	n: memoryview = memoryview(bytearray())

# FSE decoding table
class FSEDT(HDT):
	# next state
	t: memoryview = memoryview(bytearray()).cast("H")

class DZstdState:
	# byte
	b: int = 0
	# out byte
	y: int = 0
	# dictionary ID
	d: int = 0
	# window
	w: memoryview = memoryview(bytearray())
	# max block size
	m: int = 0
	# uncompressed size
	u: int = 0
	# has checksum
	c: int = 0
	# offsets
	o: memoryview = memoryview(bytearray()).cast("i")
	# window head
	e: int = 0
	# last huffman decoding table
	h: Union[None, HDT] = None
	# last FSE decoding tables
	t: Union[None, list[FSEDT]] = None
	# last block
	e: int = 0

def slc(v: memoryview, s: Union[None, int] = None, e: Union[None, int] = None):
	if (s == None or s < 0):
		s = 0
	if (e == None or e > len(v)):
		e = len(v)
	return memoryview(bytearray(bytes(v[s:e])))

def fill(v: memoryview, n: int, s: Union[None, int] = None, e: Union[None, int] = None):
	if (s == None or s < 0):
		s = 0
	if (e == None or e > len(v)):
		e = len(v)
	v[s:e] = n.to_bytes(1, byteorder="little") * (e - s)
	return v

def cpw(v: memoryview, t: int, s: Union[None, int] = None, e: Union[None, int] = None):
	if (s == None or s < 0):
		s = 0
	if (e == None or e > len(v)):
		e = len(v)
	copy_length = e - s
	v[t:t + copy_length] = v[s:s + copy_length]
	return v

ZstdErrorCode = {
	"InvalidData" : 0,
	"WindowSizeTooLarge" : 1,
	"InvalidBlockType" : 2,
	"FSEAccuracyTooHigh" : 3,
	"DistanceTooFarBack" : 4,
	"UnexpectedEOF" : 5,
}

ec = [
	'invalid zstd data',
	'window size too large (>2046MB)',
	'invalid block type',
	'FSE accuracy too high',
	'match distance too far back',
	'unexpected EOF',
]

def err(ind: int):
	raise Exception(ec[ind])

def rb(d: memoryview, b: int, n: int):
	o = 0
	for i in range(n):
		o |= d[b + i] << (i << 3)
	return o

def b4(d: memoryview, b: int):
	return (d[b] | (d[b + 1] << 8) | (d[b + 2] << 16) | (d[b + 3] << 24)) | 0

# read Zstandard frame header
def rzfh(dat: memoryview, w: Union[memoryview, int] = 1) -> Union[int, DZstdState]:
	n3 = dat[0] | (dat[1] << 8) | (dat[2] << 16)
	if n3 == 0x2FB528 and dat[3] == 253:
		# Zstandard
		flg = dat[4]
		# single segment
		ss = (flg >> 5) & 1
		# checksum
		cc = (flg >> 2) & 1
		# dict flag
		df = flg & 3
		# frame content flag
		fcf = flg >> 6
		if flg & 8 != 0:
			err(0)
		# byte
		bt = 6 - ss
		# dict bytes
		db = 4 if (df == 3) else df
		# dictionary id
		di = rb(dat, bt, db)
		bt += db
		# frame size bytes
		fsb = (1 << fcf) if (fcf != 0) else ss
		# frame source size
		fss = rb(dat, bt, fsb) + (256 if (fcf == 1) else 0)
		# window size
		ws = fss
		if ss == 0:
			wb = 1 << (10 + (dat[5] >> 3))
			ws = wb + (wb >> 3) * (dat[5] & 7);
		if ws > 2145386496:
			err(1)
		buf = memoryview(bytearray(((fss if (fss != 0) else ws) if (w == 1) else (0 if (w != 1) else ws)) + 12))
		buf[0] = 1
		buf[4] = 4
		buf[8] = 8
		s = DZstdState()
		s.b = bt + fsb
		s.y = 0
		s.l = 0
		s.d = di
		s.w = w if (w != None and w != 0 and w != 1) else buf[12:]
		s.e = ws
		s.o = buf[0:0 + (3 * 4)].cast("i")
		s.u = fss
		s.c = cc
		s.m = min(131072, ws)
		return s
	elif ((n3 >> 4) | (dat[3] << 20)) == 0x184D2A5:
		return b4(dat, 4) + 8
	err(0)

# most significant bit for nonzero
def msb(val: int):
	bits = 0
	while (1 << bits) <= val:
		bits += 1
	return bits - 1

# read finite state entropy
def rfse(dat: memoryview, bt: int, mal: int) -> tuple[int, FSEDT]:
	# table pos
	tpos = (bt << 3) + 4
	# accuracy log
	al = (dat[bt] & 15) + 5
	if (al > mal):
		err(3)
	# size
	sz = 1 << al
	# probabilities
	probs = sz
	# symbols
	sym = -1
	# repeat
	re = -1
	# high threshold
	ht = sz
	# optimization: single allocation is much faster
	buf = memoryview(bytearray(512 + (sz << 2)))
	freq = buf[0:0 + (256 * 2)].cast("h")
	# same view as freq
	dstate = buf[0:0 + (256 * 2)].cast("H")
	nstate = buf[512:512 + (sz * 2)].cast("H")
	bb1 = 512 + (sz << 1)
	syms = buf[bb1:bb1 + sz]
	nbits = buf[bb1 + sz:]
	while sym < 255 and probs > 0:
		bits = msb(probs + 1)
		cbt = tpos >> 3
		# mask
		msk = (1 << (bits + 1)) - 1
		val = 0
		if cbt < len(dat):
			val |= dat[cbt]
		if cbt + 1 < len(dat):
			val |= (dat[cbt + 1] << 8)
		if cbt + 2 < len(dat):
			val |= (dat[cbt + 2] << 16)
		val >>= (tpos & 7)
		val &= msk
		# mask (1 fewer bit)
		msk1fb = (1 << bits) - 1
		# max small value
		msv = msk - probs - 1
		# small value
		sval = val & msk1fb
		if sval < msv:
			tpos += bits
			val = sval
		else:
			tpos += bits + 1
			if val > msk1fb:
				val -= msv
		sym += 1
		val -= 1
		freq[sym] = val
		if val == -1:
			probs += val
			ht -= 1
			syms[ht] = sym
		else:
			probs -= val
		if val == 0:
			while re == 3:
				# repeat byte
				rbt = tpos >> 3
				re = ((dat[rbt] | (dat[rbt + 1] << 8)) >> (tpos & 7)) & 3
				tpos += 2
				sym += re
	if sym > 255 or probs != 0:
		err(0)
	sympos = 0
	# sym step (coprime with sz - formula from zstd source)
	sstep = (sz >> 1) + (sz >> 3) + 3
	# sym mask
	smask = sz - 1
	for s in range(sym + 1):
		sf = freq[s]
		if sf < 1:
			dstate[s] = -sf
			continue
		# This is split into two loops in zstd to avoid branching, but as JS is higher-level that is unnecessary
		for i in range(sf):
			syms[sympos] = s
			sympos = (sympos + sstep) & smask
			while sympos >= ht:
				sympos = (sympos + sstep) & smask;

	# After spreading symbols, should be zero again
	if sympos != 0:
		err(0)
	for i in range(sz):
		# next stae
		ns = dstate[syms[i]]
		dstate[syms[i]] += 1
		# num bits
		nb = nbits[i] = al - msb(ns)
		nstate[i] = (ns << nb) - sz
	FSEDT_ret = FSEDT()
	FSEDT_ret.b = al
	FSEDT_ret.s = syms
	FSEDT_ret.n = nbits
	FSEDT_ret.t = nstate
	return ((tpos + 7) >> 3, FSEDT_ret)

# read huffman
def rhu(dat: memoryview, bt: int) -> tuple[int, HDT]:
	# weight count
	wc = -1
	# buffer
	buf = memoryview(bytearray(292))
	# header byte
	hb = dat[bt]
	# huffman weights
	hw = buf[0:256]
	# rank count
	rc = buf[256:268]
	# rank index
	ri = buf[268:].cast("H")
	# NOTE: at this point bt is 1 less than expected
	if hb < 128:
		# end byte, fse decode table
		rfse_tmp = rfse(dat, bt + 1, 6)
		ebt = rfse_tmp[0]
		fdt = rfse_tmp[1]
		bt += hb
		epos = ebt << 3
		# last byte
		lb = dat[bt]
		if lb == 0:
			err(0)
		# state1
		st1 = 0
		# state2
		st2 = 0
		# state1 bits
		btr1 = fdt.b
		# state2 bits
		btr2 = btr1
		# fse pos
		# pre-increment to account for original deficit of 1
		bt += 1
		fpos = (bt << 3) - 8 + msb(lb)
		while True:
			fpos -= btr1
			if fpos < epos:
				break
			cbt = fpos >> 3
			st1 += ((dat[cbt] | (dat[cbt + 1] << 8)) >> (fpos & 7)) & ((1 << btr1) - 1)
			wc += 1
			hw[wc] = fdt.s[st1]
			fpos -= btr2
			if fpos < epos:
				break
			cbt = fpos >> 3
			st2 += ((dat[cbt] | (dat[cbt + 1] << 8)) >> (fpos & 7)) & ((1 << btr2) - 1)
			wc += 1
			hw[wc] = fdt.s[st2]
			btr1 = fdt.n[st1]
			st1 = fdt.t[st1]
			btr2 = fdt.n[st2]
			st2 = fdt.t[st2]
		wc += 1
		if wc > 255:
			err(0)
	else:
		wc = hb - 127
		# (let i = 0; i < wc; i += 2)
		for i in range(0, wc, 2):
			bt += 1
			byte = dat[bt]
			hw[i] = byte >> 4
			hw[i + 1] = byte & 15
		bt += 1
	# weight exponential sum
	wes = 0
	for i in range(wc):
		wt = hw[i];
		# bits must be at most 11, same as weight
		if wt > 11:
			err(0)
		wes += (1 << (wt - 1)) if (wt != 0) else 0
	# max bits
	mb = msb(wes) + 1
	# table size
	ts = 1 << mb
	# remaining sum
	rem = ts - wes
	# must be power of 2
	if rem & (rem - 1) != 0:
		err(0)
	hw[wc] = msb(rem) + 1
	wc += 1
	for i in range(wc):
		wt = hw[i]
		hw[i] = (mb + 1 - wt) if (wt != 0) else 0
		rc[hw[i]] += 1
	# huf buf
	hbuf = memoryview(bytearray(ts << 1))
	# symbols
	syms = hbuf[0:ts]
	# num bits
	nb = hbuf[ts:]
	ri[mb] = 0
	# (let i = mb; i > 0; i -= 1)
	for i in range(mb, 0, -1):
		pv = ri[i]
		ri[i - 1] = pv + rc[i] * (1 << (mb - i))
		fill(nb, i, pv, ri[i - 1])
	if ri[0] != ts:
		err(0)
	for i in range(wc):
		bits = hw[i]
		if bits != 0:
			code = ri[bits]
			ri[bits] = code + (1 << (mb - bits))
			fill(syms, i, code, ri[bits])
	hdt_ret = HDT()
	hdt_ret.n = nb
	hdt_ret.b = mb
	hdt_ret.s = syms
	return (bt, hdt_ret)

# Tables generated using this:
# https://gist.github.com/101arrowz/a979452d4355992cbf8f257cbffc9edd

# default literal length table
dllt = rfse(memoryview(bytearray([81, 16, 99, 140, 49, 198, 24, 99, 12, 33, 196, 24, 99, 102, 102, 134, 70, 146, 4])), 0, 6)[1]

# default match length table
dmlt = rfse(memoryview(bytearray([33, 20, 196, 24, 99, 140, 33, 132, 16, 66, 8, 33, 132, 16, 66, 8, 33, 68, 68, 68, 68, 68, 68, 68, 68, 36, 9])), 0, 6)[1]

# default offset code table
doct = rfse(memoryview(bytearray([32, 132, 16, 66, 102, 70, 68, 68, 68, 68, 36, 73, 2])), 0, 5)[1]

# bits to baseline
def b2bl(b: memoryview, s: int):
	blen = len(b)
	bl = memoryview(bytearray(blen * 4)).cast("i")
	for i in range(blen):
		bl[i] = s
		s += 1 << b[i]
	return bl

# literal length bits
llb = memoryview(bytearray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]))

# literal length baseline
llbl = b2bl(llb, 0)

# match length bits
mlb = memoryview(bytearray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]))

# match length baseline
mlbl = b2bl(mlb, 3)

# decode huffman stream
def dhu(dat: memoryview, out: memoryview, hu: HDT):
	blen = len(dat)
	ss = len(out)
	lb = dat[blen - 1]
	msk = (1 << hu.b) - 1
	eb = -hu.b
	if lb == 0:
		err(0)
	st = 0
	btr = hu.b
	pos = (blen << 3) - 8 + msb(lb) - btr
	i = -1
	while pos > eb and i < ss:
		cbt = pos >> 3
		val = 0
		if cbt < len(dat):
			val |= dat[cbt]
		if cbt + 1 < len(dat):
			val |= dat[cbt + 1] << 8
		if cbt + 2 < len(dat):
			val |= dat[cbt + 2] << 16
		val >>= pos & 7
		st = ((st << btr) | val) & msk
		i += 1
		out[i] = hu.s[st]
		btr = hu.n[st]
		pos -= btr
	if (pos != eb) or (i + 1 != ss):
		err(0)

# decode huffman stream 4x
# TODO: use workers to parallelize
def dhu4(dat: memoryview, out: memoryview, hu: HDT):
	bt = 6
	ss = len(out)
	sz1 = (ss + 3) >> 2
	sz2 = sz1 << 1
	sz3 = sz1 + sz2
	new_bt = bt
	new_bt += dat[0] | (dat[1] << 8)
	dhu(dat[bt:new_bt], out[0:sz1], hu)
	bt = new_bt
	new_bt += dat[2] | (dat[3] << 8)
	dhu(dat[bt:new_bt], out[sz1:sz2], hu)
	bt = new_bt
	new_bt += dat[4] | (dat[5] << 8)
	dhu(dat[bt:new_bt], out[sz2:sz3], hu)
	bt = new_bt
	dhu(dat[bt:], out[sz3:], hu)

# read Zstandard block
def rzb(dat: memoryview, st: DZstdState, out: Union[memoryview, None] = None):
	bt = st.b
	# byte 0
	b0 = dat[bt]
	# block type
	btype = (b0 >> 1) & 3
	st.l = b0 & 1
	sz = (b0 >> 3) | (dat[bt + 1] << 5) | (dat[bt + 2] << 13)
	# end byte for block
	bt += 3
	ebt = (bt) + sz
	if btype == 1:
		if bt >= len(dat):
			return None
		st.b = bt + 1
		if out != None:
			st.y += sz
			fill(out, dat[bt], st.y, st.y)
			return out
		return fill(memoryview(bytearray(sz)), dat[bt])
	if ebt > len(dat):
		return None
	if btype == 0:
		st.b = ebt
		if out != None:
			dat_slice = dat[bt:ebt]
			out[st.y:st.y + len(dat_slice)] = dat_slice
			st.y += sz
			return out
		return slc(dat, bt, ebt)
	if btype == 2:
		# byte 3
		b3 = dat[bt]
		# lit btype
		lbt = b3 & 3
		# size format
		sf = (b3 >> 2) & 3
		# lit src size
		lss = b3 >> 4
		# lit cmp sz
		lcs = 0
		# 4 streams
		s4 = 0
		if lbt < 2:
			if sf & 1 != 0:
				bt += 1
				lss |= (dat[bt] << 4)
				bt += 1
				lss |= ((dat[bt] << 12) if (sf & 2 != 0) else 0)
			else:
				lss = b3 >> 3
		else:
			s4 = sf
			if sf < 2:
				bt += 1
				lss |= ((dat[bt] & 63) << 4)
				lcs = (dat[bt] >> 6)
				bt += 1
				lcs |= (dat[bt] << 2)
			elif sf == 2:
				bt += 1
				lss |= (dat[bt] << 4)
				bt += 1
				lss |= ((dat[bt] & 3) << 12)
				lcs = (dat[bt] >> 2)
				bt += 1
				lcs |= (dat[bt] << 6)
			else:
				bt += 1
				lss |= (dat[bt] << 4)
				bt += 1
				lss |= ((dat[bt] & 63) << 12)
				lcs = (dat[bt] >> 6) 
				bt += 1
				lcs |= (dat[bt] << 2)
				bt += 1
				lcs |= (dat[bt] << 10)
		bt += 1
		# add literals to end - can never overlap with backreferences because unused literals always appended
		buf = out[st.y:st.y + st.m] if (out != None) else memoryview(bytearray(st.m))
		# starting point for literals
		spl = len(buf) - lss
		if lbt == 0:
			new_bt = bt + lss
			dat_slice = dat[bt:new_bt]
			buf[spl:spl + len(dat_slice)] = dat_slice
			bt = new_bt
		elif lbt == 1:
			fill(buf, dat[bt], spl)
			bt += 1
		else:
			# huffman table
			hu = st.h
			if lbt == 2:
				hud = rhu(dat, bt)
				# subtract description length
				lcs += bt - (hud[0])
				bt = hud[0]
				st.h = hud[1]
				hu = hud[1]
			elif hu == 0:
				err(0)
			new_bt = bt + lcs
			if s4 != 0:
				dhu4(dat[bt:new_bt], buf[spl:], hu)
			else:
				dhu(dat[bt:new_bt], buf[spl:], hu)
			bt = new_bt
		# num sequences
		ns = dat[bt]
		bt += 1
		if ns != 0:
			if ns == 255:
				ns = dat[bt]
				bt += 1
				ns |= (dat[bt] << 8)
				bt += 1
				ns += 0x7F00
			elif ns > 127:
				ns = ((ns - 128) << 8)
				ns |= dat[bt]
				bt += 1
			# symbol compression modes
			scm = dat[bt]
			bt += 1
			if scm & 3 != 0:
				err(0)
			dts = [dmlt, doct, dllt]
			for i in range(3):
				md = (scm >> (((3 - i) << 1) + 2)) & 3
				if md == 1:
					# rle buf
					rbuf = memoryview(bytearray([0, 0, dat[bt]]))
					bt += 1
					fsedt_tmp = FSEDT()
					fsedt_tmp.s = rbuf[2:3]
					fsedt_tmp.n = rbuf[0:1]
					fsedt_tmp.t = rbuf[0:0 + (1 * 2)].cast("H")
					fsedt_tmp.b = 0
					dts[i] = fsedt_tmp
				elif md == 2:
					# accuracy log 8 for offsets, 9 for others
					tmp_rfse = rfse(dat, bt, 9 - ((3 - i) & 1))
					bt = tmp_rfse[0]
					dts[i] = tmp_rfse[1]
				elif md == 3:
					if st.t == None:
						err(0)
					dts[i] = st.t[i]
			st.t = dts
			bmlt = dts[0]
			boct = dts[1]
			bllt = dts[2]
			lb = dat[ebt - 1]
			if lb == 0:
				err(0)
			spos = (ebt << 3) - 8 + msb(lb) - bllt.b
			cbt = spos >> 3
			oubt = 0
			lst = ((dat[cbt] | (dat[cbt + 1] << 8)) >> (spos & 7)) & ((1 << bllt.b) - 1)
			spos -= boct.b
			cbt = spos >> 3
			ost = ((dat[cbt] | (dat[cbt + 1] << 8)) >> (spos & 7)) & ((1 << boct.b) - 1)
			spos -= bmlt.b
			cbt = spos >> 3
			mst = ((dat[cbt] | (dat[cbt + 1] << 8)) >> (spos & 7)) & ((1 << bmlt.b) - 1)
			for i in range(ns):
				llc = bllt.s[lst]
				lbtr = bllt.n[lst]
				mlc = bmlt.s[mst]
				mbtr = bmlt.n[mst]
				ofc = boct.s[ost]
				obtr = boct.n[ost]

				spos -= ofc
				cbt = spos >> 3
				ofp = 1 << ofc
				off = ofp + (((dat[cbt] | (dat[cbt + 1] << 8) | (dat[cbt + 2] << 16) | (dat[cbt + 3] << 24)) >> (spos & 7)) & (ofp - 1))
				spos -= mlb[mlc]
				cbt = spos >> 3
				ml = mlbl[mlc] + (((dat[cbt] | (dat[cbt + 1] << 8) | (dat[cbt + 2] << 16)) >> (spos & 7)) & ((1 << mlb[mlc]) - 1))
				spos -= llb[llc]
				cbt = spos >> 3
				ll = llbl[llc] + (((dat[cbt] | (dat[cbt + 1] << 8) | (dat[cbt + 2] << 16)) >> (spos & 7)) & ((1 << llb[llc]) - 1))

				spos -= lbtr
				cbt = spos >> 3
				lst = bllt.t[lst] + (((dat[cbt] | (dat[cbt + 1] << 8)) >> (spos & 7)) & ((1 << lbtr) - 1))
				spos -= mbtr
				cbt = spos >> 3
				mst = bmlt.t[mst] + (((dat[cbt] | (dat[cbt + 1] << 8)) >> (spos & 7)) & ((1 << mbtr) - 1))
				spos -= obtr
				cbt = spos >> 3
				ost = boct.t[ost] + (((dat[cbt] | (dat[cbt + 1] << 8)) >> (spos & 7)) & ((1 << obtr) - 1))

				if off > 3:
					st.o[2] = st.o[1]
					st.o[1] = st.o[0]
					off -= 3
					st.o[0] = off
				else:
					idx = off - (1 if (ll != 0) else 0)
					if idx != 0:
						off = (st.o[0] - 1) if (idx == 3) else st.o[idx]
						if idx > 1:
							st.o[2] = st.o[1]
						st.o[1] = st.o[0]
						st.o[0] = off
					else:
						off = st.o[0];
				if oubt + ll > spl:
					# non-overlapping copy
					buf[oubt:oubt + ll] = buf[spl:spl + ll]
				else:
					# overlapping copy
					for i in range(ll):
						buf[oubt + i] = buf[spl + i]
				oubt += ll
				spl += ll
				stin = oubt - off
				if stin < 0:
					blen = -stin
					bs = st.e + stin
					if blen > ml:
						blen = ml
					buf[oubt:oubt + blen] = st.w[bs:bs + blen]
					oubt += blen
					ml -= blen
					stin = 0
				if oubt + ml > stin:
					# non-overlapping copy
					buf[oubt:oubt + ml] = buf[stin:stin + ml]
				else:
					# overlapping copy
					for i in range(ml):
						buf[oubt + i] = buf[stin + i]
				oubt += ml
			if oubt != spl:
				llen = len(buf) - spl
				if oubt + llen > spl:
					# non-overlapping copy
					buf[oubt:oubt + llen] = buf[spl:spl + llen]
				else:
					# overlapping copy
					for i in range(llen):
						buf[oubt + i] = buf[spl + i]
				oubt += llen
				spl += llen
			else:
				oubt = len(buf)
			if out != None:
				st.y += oubt
			else:
				buf = slc(buf, 0, oubt)
		else:
			if out != None:
				st.y += lss
				if spl != 0:
					if 0 + lss > spl:
						# non-overlapping copy
						buf[0:0 + lss] = buf[spl:spl + lss]
					else:
						# overlapping copy
						for i in range(lss):
							buf[i] = buf[spl + i]
			elif spl != 0:
				buf = slc(buf, spl)
		st.b = ebt
		return buf
	err(2)

# concat
def cct(bufs: list[memoryview], ol: int):
	if len(bufs) == 1:
		return bufs[0]
	buf = memoryview(bytearray(ol))
	b = 0
	for i in range(len(bufs)):
		chk = bufs[i]
		buf[b:b + len(chk)] = chk
		b += len(chk)
	return buf

"""
 * Decompresses Zstandard data
 * @param dat The input data
 * @param buf The output buffer. If unspecified, the function will allocate
 *            exactly enough memory to fit the decompressed data. If your
 *            data has multiple frames and you know the output size, specifying
 *            it will yield better performance.
 * @returns The decompressed data
"""
def decompress(dat: memoryview, buf: Union[None, memoryview] = None):
	bt = 0
	bufs: list[memoryview] = []
	ol = 0
	while len(dat) != 0:
		st = rzfh(dat, 1 if buf == None else buf)
		if type(st) == DZstdState:
			if buf == None:
				if len(st.w) == st.u:
					buf = st.w
					bufs.append(buf)
					ol += st.u
			else:
				bufs.append(buf)
				st.e = 0
			while st.l == 0:
				blk = rzb(dat, st, buf)
				if blk == None:
					err(5)
				if buf != None:
					st.e = st.y
				else:
					bufs.append(blk)
					ol += len(blk)
					cpw(st.w, 0, len(blk))
					w_set_pos = len(st.w) - len(blk)
					st.w[w_set_pos:w_set_pos + len(blk)] = blk
			bt = st.b + (st.c * 4)
		else:
			bt = st
		dat = dat[bt:]
	return cct(bufs, ol)

"""
 * Callback to handle data in Zstandard streams
 * @param data The data that was (de)compressed
 * @param final Whether this is the last chunk in the stream
"""
dammy = 0
# export type ZstdStreamHandler = (data: Uint8Array, final?: boolean) => unknown;

# Decompressor for Zstandard streamed data
class Decompress:
	s: Union[int, DZstdState] = -1
	c: list[memoryview] = []
	l: int = 0
	z: int = 0

	"""
	   * Creates a Zstandard decompressor
	   * @param ondata The handler for stream data
	"""
	def __init__(self, ondata = None):
		self.ondata = ondata
		self.s = None
		self.c = []
		self.l = 0
		self.z = 0

	"""
	   * Pushes data to be decompressed
	   * @param chunk The chunk of data to push
	   * @param final Whether or not this is the last chunk in the stream
	"""
	def push(self, chunk: memoryview, final: Union[bool, None] = None):
		if type(self.s) == int:
			sub = min(len(chunk), self.s)
			chunk = chunk[sub:]
			self.s -= sub
		sl = len(chunk)
		ncs = sl + self.l
		if self.s == 0 or self.s == None:
			if final is True:
				if ncs == 0:
					self.ondata(memoryview(bytearray()), True)
					return
				# min for frame + one block
				if ncs < 5:
					err(5)
			elif ncs < 18:
				self.c.append(chunk)
				self.l = ncs
				return
			if self.l != 0:
				self.c.append(chunk)
				chunk = cct(self.c, ncs)
				self.c = []
				self.l = 0
			self.s = rzfh(chunk)
			if type(self.s) == int:
				return self.push(chunk, final)
		if type(self.s) != int:
			s_zstdstate = self.s
			if ncs < (self.z if (self.z != 0) else 3):
				if final == True:
					err(5)
				self.c.append(chunk)
				self.l = ncs
				return
			if self.l != 0:
				self.c.append(chunk)
				chunk = cct(self.c, ncs)
				self.c = []
				self.l = 0
			chunk_check = 4 if (chunk[s_zstdstate.b] & 2 != 0) else (3 + ((chunk[s_zstdstate.b] >> 3) | (chunk[s_zstdstate.b + 1] << 5) | (chunk[s_zstdstate.b + 2] << 13)))
			self.z = chunk_check
			if self.z == 0 and ncs < chunk_check:
				if final == True:
					err(5)
				self.c.append(chunk)
				self.l = ncs
			else:
				self.z = 0
			while True:
				blk = rzb(chunk, s_zstdstate)
				if blk == None or len(blk) == 0:
					if final == True:
						err(5)
					adc = chunk[s_zstdstate.b:]
					s_zstdstate.b = 0
					self.c.append(adc)
					self.l += len(adc)
					return
				else:
					self.ondata(blk, False)
					cpw(s_zstdstate.w, 0, len(blk))
					w_set_pos = len(s_zstdstate.w) - len(blk)
					s_zstdstate.w[w_set_pos:w_set_pos + len(blk)] = blk
				if s_zstdstate.l != 0:
					rest = chunk[s_zstdstate.b:]
					self.s = s_zstdstate.c * 4
					self.push(rest, final)
					return
		elif final:
			err(5)

	# Handler called whenever data is decompressed
	ondata = None

# tests
if __name__ == "__main__":
	import sys
	import io
	testdata = {
		# should decode to Ok
		"Ok" : {
			"data" : b"\x28\xB5\x2F\xFD\x24\x02\x11\x00\x00\x4F\x6B\x64\x50\xA9\x5A",
			"result_data" : b"Ok",
			"result_len" : 2,
		},
		# should decode to Lorem ipsum text
		"Lorem ipsum" : {
			"data" : b"\x28\xB5\x2F\xFD\x24\x84\xFD\x02\x00\x92\x07\x15\x13\x90\x07\x0C\xC9\x6E\xBB\x5B\x93\x4C\x3E\x07\xB5\xE6\x59\xF1\x89\x6B\x52\xE9\x83\xDA\x40\xF8\x2D\xED\x1E\xA9\xCB\x01\x73\x2E\x97\x5D\xB9\x13\x5B\x66\x79\xAF\x81\x66\xE0\x43\x50\xFB\x47\xFB\x21\xFC\x89\x59\x0B\x3E\xB8\x8C\x4E\xC0\x9A\xAD\x42\x15\x72\x6D\x26\x1E\x5A\xCC\x39\xC1\x74\x72\xAE\xBD\xA6\x65\xF1\xEB\x2D\xD4\x8F\x34\x01\x01\x02\x00\x3E\x53\x53\xB3\xE6\x19\xB0\x58\x5B\x26",
			"result_data" : b"Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nullam nec sem urna. Morbi mollis, massa a convallis iaculis, mauris neque.",
			"result_len" : 132,
		},
		# should decode to 1000000 of nulls
		"1000000 null" : {
			"data" : b"\x28\xB5\x2F\xFD\xA4\x40\x42\x0F\x00\x54\x00\x00\x10\x00\x00\x01\x00\xFB\xFF\x39\xC0\x02\x02\x00\x10\x00\x02\x00\x10\x00\x02\x00\x10\x00\x02\x00\x10\x00\x02\x00\x10\x00\x02\x00\x10\x00\x03\x12\x0A\x00\xCC\xAE\xCA\x39",
			"result_len" : 1000000,
			"check_zero" : True,
		},
	}
	if "--run-tests" in sys.argv:
		for key in testdata:
			decompressed = decompress(memoryview(bytearray(testdata[key]["data"])), memoryview(bytearray(testdata[key]["result_len"])))
			result_ok = True
			if "result_data" in testdata[key] and decompressed != testdata[key]["result_data"]:
				result_ok = False
			if "result_len" in testdata[key] and len(decompressed) != testdata[key]["result_len"]:
				result_ok = False
			if "check_zero" in testdata[key] and testdata[key]["check_zero"] and sum(decompressed) != 0:
				result_ok = False
			if not result_ok:
				print(key + " test NG")
			else:
				print(key + " test OK")
	elif "--run-tests-stream" in sys.argv:
		for key in testdata:
			with io.BytesIO() as wft:
				def on_data(chunk, is_last):
					wft.write(chunk)
				decompress_obj = Decompress(on_data)
				decompress_obj.push(memoryview(bytearray(testdata[key]["data"])), False)
				decompress_obj.push(memoryview(bytearray()), True)
				wft.seek(0)
				decompressed = memoryview(bytearray(wft.read()))
				result_ok = True
				if "result_data" in testdata[key] and decompressed != testdata[key]["result_data"]:
					result_ok = False
				if "result_len" in testdata[key] and len(decompressed) != testdata[key]["result_len"]:
					result_ok = False
				if "check_zero" in testdata[key] and testdata[key]["check_zero"] and sum(decompressed) != 0:
					result_ok = False
				if not result_ok:
					print(key + " test NG")
				else:
					print(key + " test OK")
	elif len(sys.argv) > 3 and sys.argv[1] == "--test-file":
		with open(sys.argv[2], "rb") as f:
			decompressed = decompress(memoryview(bytearray(f.read())))
			with open(sys.argv[3], "wb") as wf:
				wf.write(bytes(decompressed))
	elif len(sys.argv) > 3 and sys.argv[1] == "--test-file-stream":
		with open(sys.argv[2], "rb") as f:
			with io.BytesIO() as wft:
				def on_data(chunk, is_last):
					wft.write(chunk)
				decompress_obj = Decompress(on_data)
				decompress_obj.push(memoryview(bytearray(f.read())), False)
				decompress_obj.push(memoryview(bytearray()), True)
				wft.seek(0)
				with open(sys.argv[3], "wb") as wf:
					wf.write(wft.read())
