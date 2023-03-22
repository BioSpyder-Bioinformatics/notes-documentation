# Tkinter
Tkinter is a native Python library that works on Windows/Mac/Linux to create GUIs.


## Widgets
Tkinter is based on 'widgets', eg graphical components that hold a function. The most common are:
- Label
- Button
- Entry (allows only 1 line entry)
- Text (allows multiline entry) 
- Frame (rectangle used to group related widgets)

*TKK (styled tk library) also brings in*:
- Progress bar
- Better looks
- Dropdown menu
- Much more...


## Layout managers
Tkinter has 3 main layout managers to help you lay your widgets
- Pack
	+ Pack simply shoves the widget in the first available spot. Not too nice
label.pack()

- Grid
	+ Grid places the widget in an imaginary grid, you can control the relative positioning of your widget by specifying the row and column (best imo)
label.grid(row=0, column=0 )

- Place
	+ Place places the widgets at specific x-y coordinates. 


### Simple start
Simple start of the app. You need to declare a Tk() element which is going to be 'root' or the window which pops up. From there you can start adding the widgets, then run the window in loop (so it persists on the screen)

```py
# Standard initialisation
# Create a new window
window = tk.Tk()
# Set the window title
window.title("TempoSeq Aligner")
# Set the window size
window.geometry("1200x800") #width x height

# Add a label to it
label1 = tkk.Label(window, text='Select reference genome', font=arial20)
label1.grid(row=0) # place the label(need to!!)

# Start loop 
window.mainloop() # This makes the window persist, without it the app terminates as soon as it loads all variables
```


### Using LabelFrame
When using grid, I find very useful dividing the page in quadrants which enclose sets of widgets. You can use LabelFrame to do this. Need to append the frame to the main 'window', then you can just append widgets to the labelframe object
```py
window = tk.Tk()

# Declare frames and assign them to grid!
frame1 = tk.LabelFrame(window, height=680, width=600, bg='red') 
frame1.grid(row=0, column=0)
# Other frame, these divide the page in 2
frame2 = tk.LabelFrame(window, height=680, width=600, bg='blue')
frame2.grid(row=0, column=1)

# Assign a label per frame
label1 = tkk.Label(frame1, text='Select reference genome', font=arial20)
label1.grid(row=0)

label2 = tkk.Label(frame2, text='Select reference genome', font=arial20)
label2.grid(row=0)

```
 



### Font
You can setup font objects so to have a standard styling eg
```py
# Set up font objects
arial20 = tkFont.Font(family='Arial', size=20)

# Use it on a label
label1 = tkk.Label(frame1, text='Select reference genome', font=arial20)

```




### Dropdown menu styling
Dropdown menus widgets are found in TKK and are called OptionMenu. These widgets cannot be styled directly and need to be declared, then the object has to be accessed to interact with the styling using the .config:
```py
# reference genome variable (Required to retrieve the choice option later on)
variable1 = tk.StringVar(frame1, genome_options[0]) 
# Declare variable
dropdown1 = tk.OptionMenu(frame1, variable1, *genome_options)
# Change style
dropdown1.config(font=arial15, justify='center', width=20)
# Place in grid
dropdown1.grid(row=1)
```


### Display text on Text
(Usually to input text, can be hijacked to output)
```py
# Textbox variable to be updated
text5 = tk.Text(frame2)#, width=50, height=50)
# Don't allow input!
text5['state'] = 'disabled'


# Function assigned to a button. Gets input from somewhere else and outputs the result in text box
def get_directory_or_files():
    # Connected to button, prompts the user to get the files
    var = askopenfilenames() or askdirectory() 
    # Format text to display nicely
    text = ' \n'.join([x for x in var])
    # Enable editing of text5
    text5['state'] = 'normal'
    # Set new text5 (1.0 is starting row 1, column 0; up to tk.END (end of text); text to set)
    text5.replace('1.0', tk.END, text)
    # Disable editing
    text5['state'] = 'disabled'

```




### Get input file


