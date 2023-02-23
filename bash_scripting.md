# Bash cheatsheet

### Most popular commands

- `ls`: List directory contents.
   - `-l`: Long listing format
   - `-a`: Show hidden files
   - `-h`: Show file sizes in human-readable format
   - `-t`: Sort files by modification time

- `cd`: Change the current working directory.

- `pwd`: Print the current working directory.

- `mkdir`: Create a new directory.
   - `-p`: Create parent directories if they don't exist
   * Ex: `mkdir -p /home/user/Documents/Folder1/Folder2` creates the directories `/home/user/Documents/Folder1/Folder2` if they don't exist.

- `rm`: Remove files or directories.
   - `-r`: Recursively remove a directory and its contents
   - `-f`: Force removal without confirmation
   - `-i`: Interactive mode, ask for confirmation before each removal
   * Ex: `rm -r /home/user/Documents/Folder1` removes the directory `/home/user/Documents/Folder1` and all its contents.

- `cp`: Copy files or directories.
   - `-r`: Copy directories and their contents recursively
   - `-i`: Interactive mode, ask for confirmation before overwriting an existing file
   - `-u`: Copy only when the source file is newer than the destination file
   * Ex: `cp -r /home/user/Documents/Folder1 /home/user/Documents/Folder2` copies the directory `/home/user/Documents/Folder1` to `/home/user/Documents/Folder2`.

- `mv`: Move or rename files or directories.
   - `-i`: Interactive mode, ask for confirmation before overwriting an existing file
   * Example: `mv /home/user/Documents/Folder1 /home/user/Documents/Folder2` renames the directory `/home/user/Documents/Folder1` to `/home/user/Documents/Folder2`.

- `sudo`: Execute a command with administrative privileges.
   - `-s`: Run a shell with administrative privileges

- `apt`: Manage software packages.
   - `install`: Install a package
   - `remove`: Remove a package
   - `update`: Update the package list
   - `upgrade`: Upgrade installed packages to the latest version

- `grep`: Search for a pattern in a file or output.
    - `-i`: Case-insensitive search
    - `-n`: Show line numbers of matching lines
    - `-r`: Recursively search files in subdirectories
    * Example: `grep -i "error" /var/log/syslog` searches for the word "error" in the `/var/log/syslog` file.

- `cat`: Concatenate and print files to the terminal.
    - `-n`: Number the output lines

   Example: `cat -n /etc/apt/sources.list` prints the contents of the `/etc/apt/sources.list` file with line numbers.

- `echo`: Print text to the terminal.



## Bash Scripting

*Variables*
```sh
variable_name="value"          # assign a value to a variable
echo $variable_name            # print the value of a variable
readonly variable_name="value" # create a read-only variable
unset variable_name            # delete a variable
```

_Variable types_
```sh
string="Hello, world!"         # a string variable
integer=42                     # an integer variable
float=3.14                     # a floating-point variable
array=(value1 value2 value3)   # an array variable
```


*Input/Output*
```sh
echo "Hello, world!"           # print a string to the terminal
echo -n "Enter your name: "    # print a prompt without a new line
read input_variable            # read input from the user
echo $input_variable           # print the value of the input variable
```


*Conditionals*
```sh
if [ condition ]
then
    # code to execute if the condition is true
elif [ condition ]
then
    # code to execute if the second condition is true
else
    # code to execute if both conditions are false
fi

# Conditional execution of commands
if [ $1 -gt 10 ]; then
    echo "$1 is greater than 10"
else
    echo "$1 is less than or equal to 10"
fi

```


*For-While Loops*
```sh
# FOR loops
for variable in values
do
    # code to execute for each value
done

#Ex
for (( i=0; i<10; i++ ))
do
    # code to execute for each iteration
done

#Other
for i in {1..5}; do
    echo $i
done


# WHILE loops
while [ condition ]
do
    # code to execute while the condition is true
done

# Ex
until [ condition ]
do
    # code to execute until the condition is true
done
```


*Command-Line/Function Arguments*
```sh
$0        # script name
$1        # first argument
$2        # second argument
$@        # all arguments
$#        # number of arguments
```


*File I/O*
```sh
cat file_name                    # print the contents of a file
echo "text" > file_name          # write text to a file (overwrite)
echo "text" >> file_name         # append text to a file
grep "pattern" file_name         # search for a pattern in a file
sed "s/old/new/g" file_name      # replace text in a file

# Advanced
command 2>&1                     # redirect standard error to standard output
command > /dev/null 2>&1         # discard all output and errors
command < input.txt				 # Redirecting input from a file
#ex:
sort < input.txt				 # Will sort the lines in the input file
command1 | command2				 # Piping output to another command

$(command)                       # execute a command and use its output as a value
${variable:-default_value}       # use a default value if a variable is not set
${variable:=default_value}       # set a default value if a variable is not set
${variable:+value_if_set}        # use a value if a variable is set
${variable#prefix}              # remove a prefix from a variable
${variable%suffix}              # remove a suffix from a variable



```


*Math*
```sh
result=$(expr 4 + 2)      # evaluate an expression and store the result in a variable
echo $result              # prints "6"

result=$((4 + 2))         # evaluate an arithmetic expression and store the result in a variable
echo $result              # prints "6"

result=$(bc <<< "scale=2; 1.0 / 3.0")  # evaluate a floating-point expression and store the result in a variable
echo $result                           # prints "0.33"
```


*Functions*
```sh
function_name() {
    # code to execute in the function
    return value            # optional return value
}

function_name argument    # call a function with an argument

## Usage
# define a function that takes two arguments and returns their sum
sum() {
    local result=$(($1 + $2))  # perform the arithmetic operation and store the result in a local variable
    echo $result               # print the result
}

# call the function and pass two arguments
sum 3 5  # prints "8"


# !!! Capture return value
result=$(sum 2 3)
echo $result  # prints "5"
```


_Default Argument Values_
```sh
print_greeting() {
    local name=${1:-"world"}  # if no argument is provided, use "world" as the default value
    echo "Hello, $name!"
}

print_greeting  # prints "Hello, world!"
print_greeting "Alice"  # prints "Hello, Alice!"
```

_Variable Scoping_
```sh
greeting="Hello"

print_greeting() {
    local greeting="Hi"
    echo "$greeting, $1!"
}

print_greeting "Alice"  # prints "Hi, Alice!"
echo "$greeting"  # prints "Hello"
```

_Recursion_
```sh
factorial() {
    if [ $1 -eq 1 ]; then
        echo 1
    else
        local subresult=$(factorial $(($1 - 1)))
        echo $(($1 * $subresult))
    fi
}

result=$(factorial 5)
echo $result  # prints "120"
```


*Conditional execution*
```sh
command1 && command2  # execute command2 only if command1 succeeds (returns 0)
command1 || command2  # execute command2 only if command1 fails (returns non-zero)
command1 ; command2   # execute command2 regardless of the success or failure of command1
```





### Get command line arguments
The getopts command is used in a loop to parse the command line arguments. Each time through the loop, getopts extracts the next option (if any) from the command line and sets the variable $OPTARG to the argument (if any) for that option.

```sh
#!/bin/bash

while getopts ":a:b:" opt; do
  case ${opt} in
    a )
      echo "Option -a was triggered with argument: ${OPTARG}"
      ;;
    b )
      echo "Option -b was triggered with argument: ${OPTARG}"
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." 1>&2
      exit 1
      ;;
  esac
done


# --------
$ ./script.sh -a foo -b bar
>>> Option -a was triggered with argument: foo
>>> Option -b was triggered with argument: bar


```



### Process management
Bash provides several commands for managing processes, including ps, kill, and jobs.

- `ps` displays information about running processes. By default, it shows information about the processes running in the current terminal session.
- `kill` sends a signal to a process, causing it to terminate. The most common signal is SIGTERM, which politely asks the process to terminate. If the process does not terminate, you can use SIGKILL to force it to terminate.
- `jobs` displays a list of the jobs running in the current terminal session. You can use bg and fg to move a job to the background or foreground, respectively.




### String manipulation
Bash allows to manipulate strings, the most common used commands include:
- `grep`: searches for a pattern in a file or input
	+ `grep pattern file.txt`
- `sed`: can be used to insert/delete/search/replace text in a file or input
	+ `sed 's/unix/linux/ geekfile.txt`
		* This will replace unix with linux
		* 's' stands for substitute
		* By default it only replaces the first occurrence in each line
			- `sed 's/unix/linux/1' geekfile.txt`
			- This will only replace the '1'st occurrence. Using g instead of 1 replaces all (global)
- `cut`: can be used to extract a substring from a file or input
	+ `cut -c 1-5 file.txt`
- `tr`: used to translate or delete characters
	+ `echo "Hello, world!" | tr '[:lower:]' '[:upper:]'`
		* This inverts the lower case with upper
		* Basically you pipe in input, find the first expression (eg [:lower:]) and replace with second (eg [:upper:])
- `awk`: awk is a scripting language used for manipulating data and generating reports. 
	+ It does scan files by line, split lines, compare patterns and perform actions on matches
	+ `awk options 'selection _criteria {action }' input-file > output-file`
	+ ex print only lines that have the word 'manager' in them
		* `awk '/manager/ {print}' employee.txt`
	+ ex for each line, print the first and fourth argument
		* `awk '{print $1,$4}' employee.txt`
	+ Has a lot of keyword used for operations, very useful!!





### Debugging
Bash provides several tools for debugging scripts, including the `set` command and the -x option.
- `set -e` at the top of the script causes the script to exit immediately if any command returns a non-zero exit status
- `set -x` prints each command as it is executed for monitoring
- `set +x` inverts -x





























































