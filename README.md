# psnd
Safe compress continuous functions by p0.

# Context
There's many programs to use predictor applying into compression.
So this is another try.

# Usage
    ./psnd -e < continuous.dat-16bit > compress.dat
    ./psnd -d < compress.dat > continuous.dat-16bit

# Description
Compress continuous 16-bit signed integer input by using P0 prediction. This uses block and count by factorial but this is safe for compress and uncompress.
There exists unsafe method like this without p0, it accesses aleph_\]0, 1\[ space on counting aleph_2 as infinity. It's unsafe. (cf. 1999 bit operation compression patents, 2005-2007 journals around physics but google doesn't find them that I saw.)

# Tips
This is only useful if original datastream is continuous enough, otherwise, this expands data (becase of this, this is safe.)

# Another Download Sites
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU

# Archive
This repository is archived, so without bug report, will no change.
