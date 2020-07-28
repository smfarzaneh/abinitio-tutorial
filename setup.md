## Setting up the softwares
Here is a step-by-step guide on how to install and do a test-run of [quantum ESPRESSO](https://www.quantum-espresso.org/) on a Linux machine. 

## 0. Operating system 
Any operating system (Windows, Mac OS, Linux) will do as long as a UNIX environment with a command shell is available. If you want to use Windows you are going to need *Cygwin* which mimics a UNIX environment on Windows. However, I **recommend** using **Ubuntu 18** on a virtual machine as I have seen relatively fewer issues on this OS regarding compilers and the libraries compatibility. 


## 1. Install compilers
You will need to install `C` and `Fortran` compilers. The recommended version for both is `GNU 7.5.0`.  
You can install these on Ubuntu via the following commands: 

### update the package list  

```bash
sudo apt update
```

### install `build-essentials` and `manpages-dev` packages 
```bash
sudo apt install build-essential
sudo apt install manpages-dev
```
These packages install the required compilers such as `gcc`, `g++`, `make`, `gfortran`, and etc with the corresponding version `GNU 7.5.0`. 
To make sure all is installed correctly, check the versions: 
```bash
gcc --version
gfortran --version
```

## 2. Install parallel compilers 
Next you need to install Open MPI which provides parallel libraries and compilers such as `mpif90` which will be used by Quantum ESPRESSO. Here is how to install it in Linux. 

```bash
sudo apt update
sudo apt install openmpi-bin
```
Checking the `mpif90` version, you should get the same version as that of `gfortran` that is `GNU 7.5.0`. 
```bash
mpif90 --version 
```

## 3. Download latest release of Quantum ESPRESSO
First you need to download the latest stable release of [quantum ESPRESSO](https://github.com/QEF/q-e) (as of now `6.5` version is the latest one.) 
Then you will unzip the files and go into the the directory that was just created. 
```bash
wget https://github.com/QEF/q-e/releases/download/qe-6.5/qe-6.5-ReleasePack.tgz
tar -zxvf qe-6.5-ReleasePack.tgz 
cd qe-6.5
```

## 4. Install Quantum ESPRESSO 
If all the required compilers are installed correctly and the versions are compatible, then the following two commands should run without any errors. 

### run the configure 
make sure your current directory is `qe-6.5`. 
```bash 
./configure
```

In the log that follows, the following line should show up: 
```
parallel environment detected successfully. 
```

### compile the project 
This will take a bit of time. This command will compile all the files in the project and create executables. 
```bash
make all 
```

### setting up **environment variables**
Environment variables usually describe the directories that the software works with. They are also used to configure the conditions under which the executables run. For example, you can toggle between serial and parallel processing modes. Although the environment variables are predefined, you may want to change them depending on your specific needs. 
You can find the `environment_variables` file in the root directory of quantum ESPRESSO, that is `qe-6.5/environment_variables`.
In this file, you can toggle between serial and parallel processing by commenting and uncommenting two lines of code.  
For serial processing, you leave the `PARA_PREFIX` variable empty 
```bash 
PARA_PREFIX=" "
# PARA_PREFIX="mpirun -np 4"
```
and for parallel processing, you set `PARA_PREFIX` to include a parallel prefix as follows 
```bash 
# PARA_PREFIX=" "
PARA_PREFIX="mpirun -np 4"
```
where `mpirun -np 4` means that the `mpirun` runs your executables in parallel and the number `4` denotes the *number of processors* which `-np` stands for. 

## 5. Run Quantum ESPRESSO examples 
Change directory to the `example01` folder. 
```bash 
cd PW/examples/example01
```

Run the example. 
```bash 
bash run_example
```

If everything goes smoothly and without error, at the end you should get 
```bash
done.
```