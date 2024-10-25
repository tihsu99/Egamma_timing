import subprocess, os
from termcolor import cprint
import optparse, argparse

URL_name = "root://se01.grid.nchc.org.tw//"

def set_url(url_name):
  URL_name = url_name

def listdir(dataset_path):
  result = subprocess.Popen("xrdfs " + URL_name + " ls " + dataset_path, shell=True, stdout=subprocess.PIPE)
  result = result.stdout.read()
  result = result.decode('utf-8').split('\n')
  final_result = []
  for string in result:
    if '' == string: continue
    final_result.append(string.strip())
  return final_result

def isdir(dataset_path):
  result = subprocess.Popen("xrdfs " + URL_name + " stat " + dataset_path, shell=True, stdout=subprocess.PIPE)
  result = result.stdout.read()
  result = result.decode('utf-8')
  return ('IsDir' in result)

def rmdir(datasetname):
  lists = listdir(datasetname)
  for element_ in lists:
    if isdir(element_):
      rmdir(element_)
      cprint("xrdfs " + URL_name + " rmdir " + element_, "yellow")
      os.system("xrdfs " + URL_name + " rmdir " + element_)
    else:
      cprint("xrdfs " + URL_name + " rm " + element_, "green")
      os.system("xrdfs " + URL_name + " rm " + element_)

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--target',   type=str)
  args = parser.parse_args()
  args_dict = vars(args)
  rmdir(args.target)
