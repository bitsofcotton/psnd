# psnd
Safe compress continuous functions by p0.

# Context
There's many programs to use predictor applying into compression.
So this is another try.

# Usage
    ./psnd -e < continuous.dat-16bit-signed > compress.dat
    ./psnd -d < compress.dat > continuous.dat-16bit-signed
    # To reproduce 16bit-sound monoral raw file, we can use with sox program:
    sox -e signed -L -b 16 -c 1 input.wav output.raw

# Description
Compress continuous 16-bit signed integer input by using P0 prediction. This uses block and count by factorial but this is safe for compress and uncompress.
There exists unsafe method like this without p0, it accesses aleph_\]0, 1\[ space on counting aleph_2 as infinity. It's unsafe. (cf. 1999 bit operation compression patents, 2005-2007 journals around physics but google doesn't find them that I saw.)

# Tips
This is only useful if original datastream is continuous enough, otherwise, this expands data (becase of this, this is safe.)

# Another Download Sites      
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing      
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU

# Real close
2023/03/13 try to real close. integrate all files into lieonn.hh.
2023/03/24 code clean.
2023/12/09 param change, realclose.
2024/12/02 taylor function improvement, reclose with this.
2025/04/19 merge latest lieonn.
2025/07/21 merge latest lieonn.

