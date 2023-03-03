from tkinter import *

root = Tk()
root.geometry("800x500")
root.title("Cancer Immunotherapy")
root.resizable(True,True)

menubar = Menu(root)
dmenu = Menu(menubar, tearoff=0)
dmenu.add_command(label="Import Data", command=lambda: print("Import Data pressed"))
dmenu.add_command(label="Clean Data", command=lambda: print("Clean Data pressed "))
dmenu.add_command(label="Split the Data", command= lambda: print("Split the Data pressed"))
# dmenu.add_command(label="Save as ....", command=lambda: print("Save as"))
dmenu.add_separator()
dmenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="Data Cancer", menu=dmenu)

mmenu = Menu(menubar, tearoff=0)
mmenu.add_command(label="Immune Checkpoint Blockade Therapy", command=lambda: print(""))
mmenu.add_command(label="Car T-Cell Therapy", command=lambda: print(""))
mmenu.add_command(label="Other Therapeutic Strategies", command=lambda: print(""))
menubar.add_cascade(label="   Therapy", menu=mmenu)

pmenu = Menu(menubar, tearoff=0)
pmenu.add_command(label="Predicting immunogenicity of anigens", command=lambda: print("Mechanistic level"))
pmenu.add_command(label="Predicting response to checkpoint inhibitors", command=lambda: print("Clinical level"))
pmenu.add_command(label="Improve", command=lambda: print("Improve from Predictions menu"))
menubar.add_cascade(label="   Predictions", menu=pmenu)

emenu = Menu(menubar, tearoff=0)
emenu.add_command(label="Checkpoint", command=lambda: print("Checkpoint menu"))
emenu.add_command(label="Validation", command=lambda: print("Validation menu"))
menubar.add_cascade(label="Checkpoint and Validation", menu=emenu)


amenu = Menu(menubar, tearoff=0)
amenu.add_command(label="copyright - Miftah Rangkuti @2023", command=lambda: print("copyright - Miftah Rangkuti @2023"))

menubar.add_cascade(label="About", menu=amenu)

root.config(menu=menubar )

root.mainloop()
