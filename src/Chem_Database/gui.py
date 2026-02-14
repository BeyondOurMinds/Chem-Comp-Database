import tkinter as tk
from tkinter import filedialog, messagebox, Button, Label, OptionMenu, ttk
from Chem_Database.database import Database
from Chem_Database.image import ImageHandler
from PIL import ImageTk


class app:
    
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Chem Database GUI")

        self.file_path = None
        self.database = None
        self.image_handler = None
        self.display_frame = None
        self.database_path_label = None
        self.current_index = 1
        self.current_photo = None

        file_frame = tk.Frame(self.root, bd=2, relief="groove")
        file_frame.pack(padx=10, pady=10)

        self.file_path_label = Label(file_frame, text="No file selected", border=2, relief="sunken")
        self.file_path_label.grid(row=0, column=1, padx=5, pady=10)
        file_select = Button(file_frame, text="Select SDF File", command=self.open_file)
        file_select.grid(row=0, column=0, padx=5, pady=10)
        load_sdf = Button(file_frame, text="Load SDF", command=self.load_sdf)
        load_sdf.grid(row=0, column=2, padx=5, pady=10)
        create_db = Button(file_frame, text="Create Database", command=self.create_database)
        create_db.grid(row=1, column=0, padx=5, pady=10)
        
        # testing drop down for molecule selection
        options = ["Molecule 1", "Molecule 2", "Molecule 3"]
        self.selected_option = tk.StringVar(file_frame)
        self.selected_option.set(options[0])
        dropdown = OptionMenu(file_frame, self.selected_option, *options)
        dropdown.grid(row=1, column=1, padx=5, pady=10)

        # testing collapse frame for molecule details
        self.details_frame = tk.Frame(self.root, bd=2, relief="groove")
        self.details_frame.pack(padx=10, pady=10)
        details_label = Label(self.details_frame, text="Filter Options")
        # adding checkbox filter options for molecule details
        self.filter_var1 = tk.BooleanVar()
        self.filter_var2 = tk.BooleanVar()
        self.filter_var3 = tk.BooleanVar()
        self.filter_var4 = tk.BooleanVar()
        filter1 = tk.Checkbutton(self.details_frame, text="Molecular Weight <= 500 Da", variable=self.filter_var1)
        filter2 = tk.Checkbutton(self.details_frame, text="Logp <= 5", variable=self.filter_var2)
        filter3 = tk.Checkbutton(self.details_frame, text="H-bond donors <= 5", variable=self.filter_var3)
        filter4 = tk.Checkbutton(self.details_frame, text="H-bond acceptors <= 10", variable=self.filter_var4)
        filter1.grid(row=1, column=0, sticky='w', padx=5, pady=5)
        filter2.grid(row=2, column=0, sticky='w', padx=5, pady=5)
        filter3.grid(row=3, column=0, sticky='w', padx=5, pady=5)
        filter4.grid(row=4, column=0, sticky='w', padx=5, pady=5)
        details_label.grid(row=0, column=0, padx=5, pady=5)
        collapse_button = Button(file_frame, text="Filter Options", command=self.toggle_details)
        collapse_button.grid(row=1, column=2, padx=5, pady=10)

        self.root.mainloop()
    
    def toggle_details(self):
        if self.details_frame.winfo_viewable():
            self.details_frame.pack_forget()
        else:
            self.details_frame.pack(padx=10, pady=10)

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
        
        if self.database is None:
            messagebox.showerror("Error", "Database not initialized!")
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
        filter_list = []
        if self.filter_var1.get():
            filter_list.append("Mol_Weight")
        if self.filter_var2.get():
            filter_list.append("LogP")
        if self.filter_var3.get():
            filter_list.append("H_Bond_Donors")
        if self.filter_var4.get():
            filter_list.append("H_Bond_Acceptors")
        
        try:
            self.database.create_database(filter_list=filter_list)
            db_path = self.database.db_file
            self.image_handler = ImageHandler(db_path)
            # Create and show display frame after successful database creation
            if self.display_frame is None:
                self.display_frame = tk.Frame(self.root, bd=2, relief="groove")
                self.display_frame.pack(padx=10, pady=10)

                # Initialize image display with first molecule
                self.current_index = 1
                img = self.image_handler.get_image_by_cd_id(self.current_index)
                self.current_photo = ImageTk.PhotoImage(img)
                self.img_display = Label(self.display_frame, image=self.current_photo)
                self.img_display.grid(row=1, column=1, padx=5, pady=10)
                
                self.database_path_label = Label(self.display_frame, text=db_path, border=2, relief="sunken")
                self.database_path_label.grid(row=0, column=1, padx=5, pady=10)

                self.next_img = Button(self.display_frame, text="Next Molecule", command=self.display_next)
                self.next_img.grid(row=2, column=2, padx=15, pady=10)

                self.prev_img = Button(self.display_frame, text="Previous Molecule", command=self.display_previous)
                self.prev_img.grid(row=2, column=0, padx=15, pady=10)
                
                db_label = Label(self.display_frame, text="Database Path:")
                db_label.grid(row=0, column=0, padx=5, pady=10)
            else:
                # Update label if display frame already exists
                if self.database_path_label is not None:
                    self.database_path_label.config(text=db_path)
            
            messagebox.showinfo("Success", "Database created successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to create database:\n{str(e)}")
    
    def display_next(self):
        if self.image_handler is None:
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            # Increment to next molecule
            self.current_index += 1
            img = self.image_handler.get_image_by_cd_id(self.current_index)
            self.current_photo = ImageTk.PhotoImage(img)
            self.img_display.config(image=self.current_photo)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display image:\n{str(e)}")

    def display_previous(self):
        if self.image_handler is None:
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            # Decrement to previous molecule
            self.current_index -= 1
            img = self.image_handler.get_image_by_cd_id(self.current_index)
            self.current_photo = ImageTk.PhotoImage(img)
            self.img_display.config(image=self.current_photo)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display image:\n{str(e)}")

if __name__ == "__main__":
    app()