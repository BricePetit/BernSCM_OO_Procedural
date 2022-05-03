# Reorganisation in Object Oriented form and Unified Modeling Language of a climate simulation code: Bern Simple Climate Model
This project was carried out as part of a Master's thesis in computer science at the Université Libre de Bruxelles. It is based on the BernSCM model available via this link: https://github.com/bernSCM/bernSCM. The objective of this project is to transform the BernSCM climate model, developed in a procedural paradigm, into an object-oriented form using the Unified Modelling Language (UML). We performed the transformation in C++ and Python. As we had to compare the models, we also developed a procedural way in C++ and Python. During the development of this project, we used Python 3.8 and C++20. 

## Utilisation
To run a simulation, there are multiple cases. The global first case is when we want to use the code in object oriented. Then we can choose to run it in C++ or in Python. If you use Windows, use the following commands. The name_run_file must be renamed by a name of files that are in the folder "runfiles". For instance, a valid name is run_c4mip_coupled.
```sh
cd BernSCM_OOP
cd Cpp
Cpp.exe < ../runfiles/name_run_file
```
```sh
cd BernSCM_OOP
cd Python
python main.py < ../runfiles/name_run_file
```
In the case of Linux, commands are the following:
```sh
cd BernSCM_OOP
cd Cpp
./Cpp < ../runfiles/name_run_file
```
```sh
cd BernSCM_OOP
cd Python
python3 main.py < ../runfiles/name_run_file
```
For the procedural implementation, we will do the same as before, but we change the first directory that becomes BernSCM_Procedural. We then obtain for Windows:
```sh
cd BernSCM_Procedural
cd Cpp
Cpp.exe < ../runfiles/name_run_file
```
```sh
cd BernSCM_Procedural
cd Python
python main.py < ../runfiles/name_run_file
```
And we obtain for Linux:
```sh
cd BernSCM_Procedural
cd Cpp
Cpp.exe < ../runfiles/name_run_file
```
```sh
cd BernSCM_Procedural
cd Python
python main.py < ../runfiles/name_run_file
```
Simulations will produce an output that will be store in the output folder (there is one in BernSCM_Procedural and BernSCM_OOP). A way to plot the value is to use Grace, a software on Linux. There is an example of plotting script. The name of this script is test_xmgrace.sh. It can be run by doing the following commands:
```sh
cd BernSCM_OOP
./test_xmgrace.sh
```
Or 
```sh
cd BernSCM_Procedural
./test_xmgrace.sh
```
## Requirements
Required programming languages:
- C++20
- Python 3.8

Required libraries:
- Numpy for Python
