# Welcome to a simple Partial Differential Equation Solver

## Folders in Your Computer
To run this code, you have to make sure you have the following folders in the place where you are going to clone this repository to:
    * main
    * lib
    * src
    * bin
    * soundtrack

## How to Run
You can open the Makefile to inspect what it has there, but the simple steps to run any file you want are the following:
    Open a terminal
    cd /folder_with_main_lib_src_bin_etc
    make name_of_file.exe
    ./bin/name_of_file.exe

Two important things: The file you want to run must be on the folder called 'main' and, assuming the name of the file is for example "test.cpp", when calling make, you should type only "make test.exe".

## How to Clean

The object files (.o) are eliminated automatically from the moment you have the executables (.exe) running. To erase these, just type
    make clean

This will also clean the lib which is composed of all the files you have in the 'src' folder. Don't worry about this, they are going to be all compiled again when you call "make name_of_file.exe" automatically.

This command also cleans the terminal by calling the bash file "Cleaner", so do not worry about the terminal clearing without any apparent reason.

