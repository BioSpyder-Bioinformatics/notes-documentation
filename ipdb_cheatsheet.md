# IPDB Command CheatSheet
pip install ipdb

### Base
*h(elp)*
Without argument, print the list of available commands. With a command name as argument, print help about that command.

*s(tep)*
Execute the current line, stop at the first possible occasion (either in a function that is called or in the current function).

*n(ext)*
Continue execution until the next line in the current function is reached or it returns.

*unt(il)*
Continue execution until the line with a number greater than the current one is reached or until the current frame returns.

*r(eturn)*
Continue execution until the current function returns.

*c(ont(inue))*
Continue execution, only stop when a breakpoint is encountered.

*a(rgs)*
Print the argument list of the current function.

*p expression*
Print the value of the expression.


### Advanced

*w(here)*
Print a stack trace, with the most recent frame at the bottom. An arrow indicates the “current frame”, which determines the context of most commands.

*d(own)*
Move the current frame one level down in the stack trace (to a newer frame).

*u(p)*
Move the current frame one level up in the stack trace (to an older frame).

*b(reak): [ ([filename:]lineno | function) [, condition] ]*
With a filename:line number argument, set a break there. If filename is omitted, use the current file. With a function name, set a break at the first executable line of that function. Without argument, list all breaks. Each breakpoint is assigned a number to which all the other breakpoint commands refer.

The condition argument, if present, is a string which must evaluate to true in order for the breakpoint to be honored.

*tbreak: [ ([filename:]lineno | function) [, condition] ]*
Temporary breakpoint, which is removed automatically when it is first hit. The arguments are the same as break.

*cl(ear): [bpnumber [bpnumber ...] ]*
With a space separated list of breakpoint numbers, clear those breakpoints. Without argument, clear all breaks (but first ask confirmation).

*disable bpnumber: [bpnumber ...]*
Disables the breakpoints given as a space separated list of breakpoint numbers. Disabling a breakpoint means it cannot cause the program to stop execution, but unlike clearing a breakpoint, it remains in the list of breakpoints and can be (re-)enabled.

*enable bpnumber: [bpnumber ...]*
Enables the breakpoints specified.

*ignore bpnumber count*
Sets the ignore count for the given breakpoint number. If count is omitted, the ignore count is set to 0. A breakpoint becomes active when the ignore count is zero. When non-zero, the count is decremented each time the breakpoint is reached and the breakpoint is not disabled and any associated condition evaluates to true.

*condition bpnumber condition*
condition is an expression which must evaluate to true before the breakpoint is honored. If condition is absent, any existing condition is removed; i.e., the breakpoint is made unconditional.

*run [args ...]*
Restart the debugged python program. If a string is supplied it is splitted with “shlex”, and the result is used as the new sys.argv. History, breakpoints, actions and debugger options are preserved. “restart” is an alias for “run”.

*l(ist): [first [,last]]*
List source code for the current file. Without arguments, list 11 lines around the current line or continue the previous listing. With one argument, list 11 lines starting at that line. With two arguments, list the given range; if the second argument is less than the first, it is a count.

