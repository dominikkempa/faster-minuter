Faster-Minuter index
====================


Description
-----------

This repository contains the implementation of a wavelet tree data
structure using the fixed block boosting (FBB) technique. The result
of plugging this wavelet tree into an FM-index is called the Faster
Minuter index.

For the description of the new wavelet tree, the Faster Minuter
index, and reports of experimental evaluation, refer to the following
paper.

    @article{fasterminuter,
      author    = {Simon Gog and Juha K{\"{a}}rkk{\"{a}}inen and
                   Dominik Kempa and Matthias Petri and Simon J. Puglisi},
      title     = {Fixed Block Compression Boosting in FM-Indexes:
                   Theory and Practice},
      journal   = {Algorithmica},
      volume    = {81},
      number    = {4},
      pages     = {1370--1391},
      year      = {2019},
      doi       = {10.1007/s00453-018-0475-9},
    }

The latest version of the wavelet tree based on FBB is available from
https://github.com/dominikkempa/faster-minuter.



Usage
-----

The wavelet tree is implemented as a C++ class wt_fbb. The class is
compatible with the sdsl library (https://github.com/simongog/sdsl-lite).
The current version has been tested on Linux/PC.

The class comes with default parameters chosen for good overall
performance and in nearly all practical scenarios should be used as
is. The key customization of the class is in plugging in different
bitvector implementations. The default is the hybrid bitvector, but
faster (and larger) alternatives are available in the sdsl library.
Refer to the paper above for details and experimental comparisons.



Terms of use
------------

wt_fbb is released under the MIT/X11 license. See the file LICENCE for
more details. If you use this code, please cite the paper mentioned
above.



Authors
-------

wt_fbb was implemented by:
- [Dominik Kempa](https://scholar.google.com/citations?user=r0Kn9IUAAAAJ)
- [Juha Karkkainen](https://scholar.google.com/citations?user=oZepo1cAAAAJ)
