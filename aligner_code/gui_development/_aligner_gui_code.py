import tkinter as tk
import tkinter.ttk as tkk # additional widgets

# Available widgets (in tk)
# Label
# Button
# Entry (allows only 1 line entry)
# Text (allows multiline entry) 
# Frame (rectangel used to group related widgets)#

# TKK also brings in
# Progress bar
# Better looks
# Much more 
# 




# Create a window -> instance of Tkinter class
window = tk.Tk()

# Add a label widget
greeting = tk.Label(text='Hello, Tkiner')

#Pass widget to window
greeting.pack()


# Label tkk
label = tkk.Label(text='HOLA', 
    foreground='blue', 
    background='red',
    width=100,
    )
label.pack()



# Make a button
button = tk.Button(
    text='Click me tk!'
)
button.pack()

button2 = tkk.Button(
    text='Tkk click me!'
)
button2.pack() #better



# Single line entry widget
entry = tk.Entry(fg='yellow', bg='blue', width=50)
entry.pack()

# Retrieve stuff with .get(), delete with .delete() and insert with .insert()


# For multi line everything is the same except for .get() that needs to get the slice of characters to retrieve (usually to get all line, in this case first line) text_box.get('1.0', tk.END) tk.END gets everything until the end
text_box = tk.Text()
text_box.pack()


# For frames you need to declare the components inside the frame directly inside it



# Send the window in loop so that it does persist
window.mainloop()

