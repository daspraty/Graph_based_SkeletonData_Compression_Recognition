
## Introduction
This repository holds the codebase, dataset and models for the paper:

**Das, Pratyusha, and Antonio Ortega. “Graph-Based Skeleton Data Compression.” 2020 IEEE
22nd International Workshop on Multimedia Signal Processing (MMSP), 21 Sept. 2020, (https://ieeexplore.ieee.org/document/9287103)


The code is written for any kinect v2 (for example NTU-RGB120 : https://github.com/shahroudy/NTURGB-D) or openpose dataset.



## Prerequisites
- Matlab R2020b
https://www.mathworks.com/products/matlab.html





### Installation
``` shell
git clone https://github.com/daspraty/skel_compression_matlab.git; 
```



please refer to below guidances to run the code for Kinect

1. run main_gft_dct_kinect.m for compressing kinect skeleton data
``` matlab main_gft_dct_kinect.m ``` 
2. Please change the source path to your datafolder.
3. Please change the output path to save your files

please refer to below guidances to run the code for openpose

1. run main_gft_dct_openpose.m for compressing openpose skeleton data
``` matlab main_gft_dct_openpose.m ``` 
2. Please change the source path to your datafolder.
3. Please change the output path to save your files

#Recognition
To run activity recognition using spatio-temporal graph based features 

1. Please run recognition.m
``` matlab recognition.m ``` 
2. Change the source data folder to the original data or compressed data folder to perform and compare recognition on the original data and compressed data
## Citation
Please cite the following paper if you use this repository in your research.
```
@INPROCEEDINGS{9287103,
  author={Das, Pratyusha and Ortega, Antonio},
  booktitle={2020 IEEE 22nd International Workshop on Multimedia Signal Processing (MMSP)}, 
  title={Graph-based skeleton data compression}, 
  year={2020},
  volume={},
  number={},
  pages={1-6},
  doi={10.1109/MMSP48831.2020.9287103}}
```

## Contact
For any question, feel free to contact
```
Pratyusha Das    : daspraty@usc.edu

```
