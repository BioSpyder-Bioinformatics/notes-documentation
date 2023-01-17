# Regular Expressions in Python
Regular Expression (RE) is a sequence of characters that forms a search pattern

Python builtin package:
`import re`

- You can search both the unicode strings (str) as well as 8-bit strings (bytes). However these cannot be mixed in patterns or in the find/replace function.
- To avoid issues with backslashes, feed in the 'raw' strings (add r before string eg r"raw string")

*For more complex functionalities* have a look a the third party package 'regex'

*TO CHECK EXPRESSION MATCHING* https://github.com/python/cpython/blob/3.11/Tools/demo/redemo.py




### *re* module commonly used methods:
`x = re.fn(r"regExpr", txt_to_be_searched)` -> Best using a raw string for the RE search as escape characters might be interpreted wrongly(r"string")
- findall
	+ Returns a _list_ containing all the matches
	+ finditer method is the same but returns an iterator     
- search
	+ Returns first occurrence of match in a _Match object_ if any (else None)
	+ match method works similar but only if at the beginning of string (useful only if applying multiline conditions)
- split
	+ Returns a _list_ where the string has been split at each match. Can add a third positional parameter to indicate the maximum number of splits
- sub
	+ Replaces one or many matches with a string. (re.sub("toFind", "toReplace", txt_to_be_searched))
	+ subn method is same as sub but returns a tuple (new_string, number_of_subs_made)
- compile('regExpr')
	+ Saves inserted regular expression for reuse, returns _Pattern object_ which is more efficient computationally
- escape(pattern)
	+ Returns regular expression escaped (just adds the \ to make sure you search exactly the inserted pattern)




### *re* module data structures
_Match object_
Match objects get returned as result of a successful search. It has several properties and methods. The most commons are:
- span() 
	+ Returns a tuple containing the start and end positions of the match
- string (property)
	+ Returns the string passed into the function (the one to be searched)
- group()
	+ Returns the part of the string where there was a match
- start() - end() - span()
	+ Returns starting/ending/both(in tuple) positions of the match



_Pattern Object_
Pattern objects are compiled RE, which have methods for various operation, including pattern matches or string replacement. 
Declaration:
`p = re.compile('ab*')`
It also accepts optional _flags_ argument, used to enable special features and syntax variations (ex re.compile('ab*', re.IGNORECASE)).
Commonly used flags:
- ASCII, A
	+ Makes escaped characters (\n, \t, ..) match only on ASCII characters with the respective properties (newline, tab, ..)
- DOTALL, S
	+ Make . match any character, including newlines
- IGNORECASE, I
	+ Do case-insensitive matches
- LOCALE, L
	+ Do a locale-aware match (??)
- MULTILINE, M
	+ Multi-line matching, affecting ^ and $
- VERBOSE, X
	+ Enable verbose RE, to organise it more cleanly (look at documentation, it's a topic on its own)

Most commonly used methods for the pattern object include:
- match()
	+ Determine if RE matches at the beginning of the string. Returns Match object if found, None otherwise
- search()
	+ Determine if RE matches at any point of the string. Returns Match object if found, None otherwise
- findall()
	+ Find all substrings where the RE matches, and returns them as a list
- finditer()
	+ Find all substrings where the RE matches, and returns them as an iterator
- split(string, maxsplit)
	+ Identical to re.split()
- sub/subn(replacement_string, string, count)
	+ Identical to re.sub/subn()
- groups
	+ Number of capturing groups in the pattern
- pattern
	+ Pattern that constitutes the object





# Components of RE strings

### Metacharacters
- \ 
	+ Signals a special sequence or the escaping of a character ("\d")
- . 
	+ Any character (except newline character) ("he..o")
- ^ 
	+ Starts with ("^hello")
- $
	+ Ends with ("hello$")
- `*`
	+ Zero or more occurrences ("he`*`o") (refers to character preceding, so this would match "ho", "heeeeeeeo", but not "he" (bc it's missing the o))
- +
	+ One or more occurrences ("he+o") (same as above)
	+ `*+`, ++, ?+
	+ The `*, +, ?` quantifiers are all greedy, as they match as much text as possible. Adding _+_ at the end of the expression will cause an expression to become *possessive* and doesn't allow backtracking.
- ?
	+ Zero or one occurrence ("he?o") (same as above)
	+ `*?`, +?, ??
		* Adding _?_ after the quantifier makes it perform in *non-greedy (or minimal)* fashion; as few characters as possible will be matched
- {n}
	+ Exactly the specified number of occurrences ("hel{2}o") (same as above, it will match "hello", but not "helo")
	+ {m,n} - {m, }
		* You can also give a range of expected occurrences, or leave the range open to mean minimum 'm', no maximum 
	+ {m,n}?
		* Causes the resulting RE to match from m to n repetitions, attempting to match as _few_ as possible
	+ {m,n}+
		* Causes the resulting RE to match from m to n repetitions, attempting to match as _many_ as possible *without backtracking*
- | 
	+ Either or, it will match one of the two specified RE ("falls|stays") (it will match exactly falls OR stay)
- ()
	+ Capture and _GROUP_, used to group character as a single unit. Useful when searching >1 pattern per expression. 
	+ For example we want to search for uppercases and numbers in the same run, we could search for uppercase "[\bA-Z+\b]" and numbers "[\b0-9+]" (\b to indicate at beginning or end of word). Combined it comes up to: `"(\b[A-Z]+\b).+(\b\d+)"`
		* The .+ section is pretty much saying that before the group we have a bunch of characters (+ indicates that) that we can ignore
		* When using () with multiple groups of characters, we can use the _Match object_ method _.group(n)_ to indicate with n which group's results we want to extract. Ex our RE on the string "The price of PINEAPPLE ice cream is 20"; result.group(1) -> 'PINEAPPLE'; result.group(2) -> 20.
	+ Making a pattern such as (?...) allows you to match the strings, but not retrieve them. I am sure there is an use for this, but it's not very obvious to me, please refer to documentation
- []
	+ A _SET_ of characters ("[a-m]"), different applications include
	+ [arn]
		* Match where one of the specified characters is present
	+ [a-n]
		* Match where there is a lower case, alphabetically between a and n
	+ [^arn]
		* Match everything EXCEPT a,r,n
	+ [0123]
		* Match where any of the specified digits are present
	+ [0-9]
		* Match any digit between 0 and 9
	+ `[0-5][0-9]`
		* Match any two digit number from 00 to 59
	+ [a-zA-Z]
		* Match any character upper or lower case
	+ [+]
		* Inside sets, `+, *, ., |, (), $, {}` have no special meaning, so [+] means return any match for + character in the string






### Special Sequences
- \A
	+ Returns a match if the following characters are at the beginning of the string ("\AThe")
- \b
	+ Returns a match where the specified characters are at the beginning or at the end of a word ("\bain" matches aint, "ain\b" matches pain)
- \B
	+ Returns a match where the characters are present *NOT* at the beginning (or end) of a word
- \d
	+ Returns a match where the string contains digits
- \D
	+ Returns a match where the string *DOESNT* contain digits
- \s
	+ Returns a match where the string contains a white space character
- \S
	+ Match where the string does *NOT* contain a white space character
- \w
	+ Returns a match where the string contains any word characters (eg a-z A-Z 0-9 and underscore character `_`)
- \W
	+ Returns a match where the string does *NOT* contain any word characters
- \Z
	+ Returns a match is the specified characters are at the end of the string




















































































