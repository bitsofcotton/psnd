# psnd
Safe compress continuous functions by p0.

# Usage
    cat continuous.dat-16-bit | ./psnd -e > compress.dat
    cat compress.dat | ./psnd -d > continuous.dat-16-bit

# Description
Compress continuous 16-bit signed integer input by using P0 prediction. This uses block and count by factorial but this is safe for compress and uncompress.
There exists unsafe method like this without p0, it accesses aleph_\]0, 1\[ space on counting aleph_2 as infinity. It's unsafe.

# Tips
This is only useful if original datastream is continuous enough, otherwise, this expands data (becase of this, this is safe.)
