# Welcome to a simple Partial Differential Equation Solver

## Folders in Your Computer
To run this code, you have to make sure you have the following folders in the place where you are going to clone this repository to:
*  _main_
*  _lib_
*  _src_
*  _bin_
*  _soundtrack_
*  _python_
*  _graphics_

## How to Run
You can open the Makefile to inspect what it has there, but the simple steps to run any file you want are the following:
* Open a terminal
* ```cd /folder_with_main_lib_src_bin_etc ```
* ```make name_of_file.exe```
* ```./bin/name_of_file.exe```

Two important things: The file you want to run must be on the folder called _main_ and, assuming the name of the file is for example __test.cpp__, when calling make, you should type only __make test.exe__.

## How to Clean

The object files (_.o_) are eliminated automatically from the moment you have the executables (_.exe_) running. To erase these, just type
* ```make clean```

This will also clean the lib which is composed of all the files you have in the _src_ folder. Don't worry about this, they are going to be all compiled again when you call "make name_of_file.exe" automatically.

This command also cleans the terminal by calling the bash file _Cleaner_, so do not worry about the terminal clearing without any apparent reason.

## Where to See the Results

If you want to see the results of your code, everything should be in the _graphics_ folder. It is there that the program __main.cpp__ will send its results, for example. This way, things will be tidier. When you run the command
* ```make py```
it will create graphics of the data given in the __output.txt__ file, this one inside the _graphics_ folder as well.

## FAQ

In the python file __graphicmaker.py__, the _filename_ variable points to the place where the __output.txt__ file should be. You must adapt this to the folder where you are going to have this repository.

