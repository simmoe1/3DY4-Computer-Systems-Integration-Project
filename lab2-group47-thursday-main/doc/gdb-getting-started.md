# GDB debugger

Debugging tools help you figure out why/where your program crashed or why the output is incorrect. They will save you time and keep you from having to write out dozens of print statements. Nevertheless, it is important to note that, while debugging tools are an essential productivity boost, they are no substitute for code reviews that help identify problems in both the logic and the implementation of your code.

GDB is a debugger that works in the terminal mode and is available by default on most Unix distributions, including the Linux version running on Raspberry Pi hardware.

# Basic usage

Assuming the source file is `my_src.cpp` and the name of the executable is `my_sw`, when compiling the code without debugging use:

`g++ my_src.cpp -o my_sw`

Running the software from the terminal can be done as:

`./my_sw`

Before using GDB add the `-g` option so that debugging info is added to the executable, for example:

`g++ -g my_src.cpp -o my_sw`

Without diving deeper, one of the most useful features of GDB is to simply find what line of code your program crashed on. You can achieve this by running your program in GDB as follows:

`gdb ./my_sw`

Then enter the following for running your program in a debug environment:

`run`

Enter the following to print the stack trace which lists all the function calls that lead to the crash. In particular look for what line(s) of code in **your files** the program crashed at after typing:

`backtrace`

When finished using GDB, exit by entering `quit`.

As shown in the references given below, GDB has more useful features, such as being able to step through your code line-by-line, printing variable values during or after your program crashes (ex: print your iterator variable in a loop to see what iteration the program crashed on), and more.

# References

If you are new to GDB and have problems getting started with it, it is suggested that you watch the first 4 mins of this [third-party video tutorial](https://www.youtube.com/watch?v=bWH-nL7v5F4) to see how to step through your code in GDB.

For another short overview to get you started, you are referred to this [textual tutorial](http://www.cs.toronto.edu/~krueger/csc209h/tut/gdb_tutorial.html).

Once more comfortable with GDB, you are referred to this short list of [gdb commands](https://www.tutorialspoint.com/gnu_debugger/gdb_commands.htm) and its
[official documentation](https://www.gnu.org/software/gdb/documentation/)
for any further clarifications.
