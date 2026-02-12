from tkinter import *
from tkinter import filedialog, messagebox
from Chem_Database.database import Database

class app:
    
    def __init__(self):
        self.root = Tk()
        self.root.title("Chem Database GUI")

        self.file_path = None
        self.database = None

        self.file_path_label = Label(self.root, text="No file selected", border=2, relief="sunken")
        self.file_path_label.grid(row=0, column=1, padx=5, pady=10)
        file_select = Button(self.root, text="Select SDF File", command=self.open_file)
        file_select.grid(row=0, column=0, padx=5, pady=10)
        load_sdf = Button(self.root, text="Load SDF", command=self.load_sdf)
        load_sdf.grid(row=0, column=2, padx=5, pady=10)
        create_db = Button(self.root, text="Create Database", command=self.create_database)
        create_db.grid(row=1, column=0, padx=5, pady=10)

        self.root.mainloop()

    def open_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("SDF files", "*.sdf")])
        if file_path:
            self.file_path = file_path
            self.file_path_label.config(text=file_path)
            # Create a new Database instance with the selected file
            self.database = Database(self.file_path)
            messagebox.showinfo("File Selected", f"Selected: {file_path}")

    def load_sdf(self):
        if self.file_path is None:
            messagebox.showerror("Error", "Please select an SDF file first!")
            return
        
        try:
            self.database.load_sdf()
            messagebox.showinfo("Success", "SDF file loaded successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load SDF file:\n{str(e)}")

    def create_database(self):
        if self.file_path is None:
            messagebox.showerror("Error", "Please select an SDF file first!")
            return
        
        if self.database is None or self.database.df is None:
            messagebox.showerror("Error", "Please load the SDF file first!")
            return
        
        try:
            self.database.create_database()
            messagebox.showinfo("Success", "Database created successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to create database:\n{str(e)}")

if __name__ == "__main__":
    app()