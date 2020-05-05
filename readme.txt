%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read-me for cl_realsimu7e12.c simulation package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Chung-Chuan Lo (cclo@mx.nthu.edu.tw)
%                 Cclolab  (http://life.nthu.edu.tw/~lablcc/index.html) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cl_realsimu7e12.c is a public version of the C code that was used 
in Kao Kuo Wei and Chung-Chuan Lo "Short Term Depression, 
Presynaptic Inhibition  and Local Neuron Diversity Play Key Functional 
Roles in the Insect Antennal Lobe". The early version of
the code was developed by Stefano Fusi and was further modified by 
Chung-Chuan Lo. If you publish works that make use of the code, 
proper acknowledgments should be given. 

The code needs to be compiled with rando2.h which is also included. 
When you run the program, it reads network.conf for the definition 
of the network structure and network.pro for the task protocol. You can 
find a list of command line arguments that you can use with the program 
in the main() section of the code. The program exports files that contain 
firing rate of each population. Depending on the arguments you use, 
the code can also export files containing spike timing, membrane potential
and conductance etc. The code was developed for general purposes and the
code provides many figures that are not used by the supplied network.conf 
and network.pro.

The included network structure file (network.conf) and task protocol 
file (network.pro) are in plain text format and are self-explained. 
The network.conf file shows an example of the artificial network applied in 
our paper "Short Term Depression, Presynaptic Inhibition  and Local Neuron 
Diversity Play Key Functional Roles in the Insect Antennal Lobe". 

Note that if any other people wants to obtain this code, please always 
forward their requests to me (Kao Kuo Wei email: kkshxt@lolab-nthu.org)
or to the lab PI Chung-Chuan Lo (cclo@mx.nthu.edu.tw).
We need to keep a track of the distribution of our codes.

How to compile the code:
We reccomand using gcc compiler
In its most basic form, gcc compiler can be used as:

gcc cl_realsimu7e12.c 

The above command executes the complete compilation process and outputs an executable with name a.out.
Use option -o, as shown below, to specify the output file name for the executable.

gcc cl_realsimu7e12.c -o cl_realsimu7e12

#############################
Kao Kuo Wei 
National Tsing Hua University 
Tel: (886) 928272399
Email: kkshxt@lolab-nthu.org
Lab website: http://life.nthu.edu.tw/~lablcc/index.html


