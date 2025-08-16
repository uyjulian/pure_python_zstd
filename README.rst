Pure Python Zstandard decoder
=============================

Pure Python Zstandard decoder. Largely based on
`fzstd <https://github.com/101arrowz/fzstd>`__.

At the moment, the decoder only works on trivial compressed data.

| **Avoid using this in production.**
| Consider using
  `zstandard <https://pypi.org/project/zstandard/>`__,
  `zstd <https://pypi.org/project/zstd/>`__, or
  `pyzstd <https://pypi.org/project/pyzstd/>`__, which has much
  better performance than this.

License
=======

MIT license. See ``LICENSE`` for more details.
