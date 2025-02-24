# Quality Assurance analysis for PET/CT and PET/MR systems with the 3D-Printed Phantom: Main Script.

import tkinter as tk
import subprocess


def select(option):
    window.destroy()
    subprocess.call(["python", option])


def execute_all():
    window.destroy()
    subprocess.call(["python", "Quantification.py"])
    subprocess.call(["python", "Resolution.py"])
    subprocess.call(["python", "Coregister.py"])
    subprocess.call(["python", "Distortion_MR.py"])
    subprocess.call(["python", "Radiomics.py"])


def execute_selected():
    selected = variable.get() + ".py"
    if selected == "Execute All":
        execute_all()
    else:
        window.destroy()
        subprocess.call(["python", selected])


window = tk.Tk()
window.title('Quality Assurance')
window.geometry('230x120')

lbl1 = tk.Label(window, text="Select type of analysis:")
lbl1.place(x=58, y=10)

options = ["Quantification", "Resolution", "Coregister", "Distortion", "Radiomics", "Execute All"]
variable = tk.StringVar(window)
variable.set(options[0])

dropdown = tk.OptionMenu(window, variable, *options)
dropdown.place(x=65, y=35)

run_button = tk.Button(window, text="Run", command=execute_selected)
run_button.place(x=100, y=70)

window.mainloop()
