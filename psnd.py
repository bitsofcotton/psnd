import subprocess
import sys

if(sys.argv[1][0] == 'e'):
  cmd = ['ksh', '-c', 'cat ' + sys.argv[3]]
  for s in range(0, int(sys.argv[2])):
    cmd[2] += ' | ./psnd -e | gzip - '
  subprocess.Popen(cmd)
elif(sys.argv[1][0] == 'd'):
  cmd = ['ksh', '-c', 'cat ' + sys.argv[3]]
  for s in range(0, int(sys.argv[2])):
    cmd[2] += ' | gunzip - | ./psnd -d '
  subprocess.Popen(cmd)
else:
  print "no such command."

