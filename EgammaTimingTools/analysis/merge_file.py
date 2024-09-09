import os, sys
import optparse, argparse


if __name__ == '__main__':
    usage = 'usage: %prog [options]'
    parser = argparse.ArgumentParser(description = usage)
    parser.add_argument('--indir', type=str)
    parser.add_argument('--particle', type=str)
    args = parser.parse_args()
    
    region = ['endcap']
    for region_ in region:
      INDIR = os.path.join(args.indir, args.particle, region_)
      samples = os.listdir(INDIR)
      for sample_ in samples:
        SAMPLE_DIR = os.path.join(INDIR, sample_)
        print('Merge {} in {}'.format(sample_, SAMPLE_DIR))
        merge_file_name = os.path.join(SAMPLE_DIR, '{}.root'.format(sample_))
        os.system('rm {}'.format(merge_file_name))
        os.system('hadd {} {}/ntupleTree_*.root'.format(merge_file_name, SAMPLE_DIR))
