import tkinter as tk
import os
from tkinter import filedialog, messagebox, Button, Label, OptionMenu
from Chem_Database.database import Database
from Chem_Database.image import ImageHandler
from Chem_Database.img_display import InfoHandler
from PIL import ImageTk, Image
import sys
import sqlite3

def resource_path(relative_path):
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)



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
        self.sort_column = "CdId"
        self.sort_direction = "ASC"

        file_frame = tk.Frame(self.root, bd=2, relief="groove")
        file_frame.grid(row=0, column=0, padx=10, pady=10)

        self.file_path_label = Label(file_frame, text="No file selected", border=2, relief="sunken")
        self.file_path_label.grid(row=0, column=1, padx=5, pady=10)

        file_select = Button(file_frame, text="Select SDF File", command=self.open_file)
        file_select.grid(row=0, column=0, padx=5, pady=10)

        load_sdf = Button(file_frame, text="Load SDF", command=self.load_sdf)
        load_sdf.grid(row=0, column=2, padx=5, pady=10)

        create_db = Button(file_frame, text="Create Database", command=self.create_database)
        create_db.grid(row=1, column=0, padx=5, pady=10)

        load_db = Button(file_frame, text="Load Database", command=self.load_database)
        load_db.grid(row=1, column=2, padx=5, pady=10)
        

        # testing collapse frame for molecule details
        self.details_frame = tk.Frame(file_frame, bd=2, relief="groove")
        self.details_frame.grid(row=3, column=0, columnspan=4, padx=10, pady=10)
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
        collapse_button.grid(row=1, column=1, padx=5, pady=10)

        self.root.mainloop()
    
    def toggle_details(self):
        if self.details_frame.winfo_viewable():
            self.details_frame.grid_remove()
        else:
            self.details_frame.grid(row=3, column=0, columnspan=4, padx=10, pady=10)

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
                self.display_frame.grid(row=1, column=0, padx=10, pady=10)


                # Initialize image display with first molecule
                # Initialize sorting defaults
                self.sort_column = "CdId"
                self.sort_direction = "ASC"
                self.current_index = 0

                # Create image label (empty initially)
                self.img_display = Label(self.display_frame)
                self.img_display.grid(row=1, column=0, padx=5, pady=10)

                self.info_display = tk.Frame(self.display_frame, bd=2, relief="groove")
                self.info_display.grid(row=1, column=1, padx=5, pady=10)

                # Load first molecule
                self.refresh_display()

                
                self.database_path_label = Label(self.display_frame, text=db_path, border=2, relief="sunken")
                self.database_path_label.grid(row=0, column=1, padx=5, pady=10)

                self.next_img = Button(self.display_frame, text="Next Molecule", command=self.display_next)
                self.next_img.grid(row=2, column=2, padx=15, pady=10)

                save_icon_path = resource_path("images/save-icon.png")
                img2 = Image.open(save_icon_path).resize((30, 30))
                self.save_icon = ImageTk.PhotoImage(img2)
                self.save_image = Button(self.display_frame, image=self.save_icon, command=self.save_current_image, width=30, height=30)
                self.save_image.grid(row=2, column=1, padx=15, pady=10)

                self.display_jump_label = Label(self.display_frame, text="Go to molecule #:", width=15)
                self.display_jump_label.grid(row=2, column=0, padx=5, pady=10)
                self.display_jump_entry = tk.Entry(self.display_frame, width=10)
                self.display_jump_entry.grid(row=2, column=4, padx=5, pady=10)

                self.chg_ord = Button(self.display_frame, text="Change Sort Order", command=self.update_order)
                self.chg_ord.grid(row=2, column=3, padx=15, pady=10)

                self.prev_img = Button(self.display_frame, text="Previous Molecule", command=self.display_previous)
                self.prev_img.grid(row=2, column=0, padx=15, pady=10)
                
                db_label = Label(self.display_frame, text="Database Path:")
                db_label.grid(row=0, column=0, padx=5, pady=10)

                options = ["CdId","Molecular weight", "LogP", "H-bond Donors", "H-bond Acceptors"]
                self.selected_option = tk.StringVar(self.display_frame)
                self.selected_option.set(options[0])
                dropdown = OptionMenu(self.display_frame, self.selected_option, *options, command=self.update_sort)
                dropdown.grid(row=0, column=2, padx=5, pady=10)
            else:
                # Update label if display frame already exists
                if self.database_path_label is not None:
                    self.database_path_label.config(text=db_path)
            
            messagebox.showinfo("Success", "Database created successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to create database:\n{str(e)}")
    
    def update_sort(self, selected_option):
        column_mapping = {
            "CdId": "CdId",
            "Molecular weight": "Mol_Weight",
            "LogP": "LogP",
            "H-bond Donors": "H_Bond_Donors",
            "H-bond Acceptors": "H_Bond_Acceptors"
        }
        self.sort_column = column_mapping.get(selected_option, "CdId")  # Default to CdId if option not found
        self.current_direction = "ASC"  # Reset to first molecule when sort order changes
        self.current_index = 0
        self.refresh_display()
    
    def update_order(self):
        # Toggle sort direction
        self.sort_direction = "DESC" if self.sort_direction == "ASC" else "ASC"
        self.current_index = 0  # Reset to first molecule when sort order changes
        self.refresh_display()
    
    def display_next(self):
        if self.image_handler is None:
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            # Increment to next molecule
            self.current_index += 1
            self.refresh_display()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display image:\n{str(e)}")

    def display_previous(self):
        if self.image_handler is None:
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            # Decrement to previous molecule
            if self.current_index > 0:
                self.current_index -= 1
                self.refresh_display()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display image:\n{str(e)}")
    
    def refresh_display(self):
        if self.image_handler is None:
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        img = self.image_handler.get_image_by_offset(
            self.current_index,
            sort_column=self.sort_column,
            sort_direction=self.sort_direction
        )

        self.current_photo = ImageTk.PhotoImage(img)
        self.img_display.config(image=self.current_photo)

        self.update_info_display()
    
    def update_info_display(self):
        # clear previous info
        for widget in self.info_display.winfo_children():
            widget.destroy()
        
        if self.database is None or self.database.db_file is None:
            messagebox.showerror("Error", "Database not initialized!")
            return
        
        info_handler = InfoHandler(self.database.db_file)
        try:
            info = info_handler.get_info_by_offset(
                self.current_index,
                sort_column=self.sort_column,
                sort_direction=self.sort_direction
            )
            for key, value in info.items():
                label = Label(self.info_display, text=f"{key}: {value}", width=60, wraplength=400, justify="left", anchor="w")
                label.pack(anchor='w')
        except Exception as e:
            messagebox.showerror("Error", f"Failed to retrieve information:\n{str(e)}")
    
    def save_current_image(self):
        if self.image_handler is None:
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            img = self.image_handler.get_image_by_offset(self.current_index)
            save_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
            if save_path:
                img.save(save_path)
                messagebox.showinfo("Success", f"Image saved to {save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save image:\n{str(e)}")
        
    def load_database(self):
        db_path = filedialog.askopenfilename(filetypes=[("SQLite Database", "*.db")])
        if db_path:
            self.database = Database(sdf_file=None, db_file=db_path)
            self.image_handler = ImageHandler(db_path)
            messagebox.showinfo("Database Loaded", f"Loaded database from: {db_path}")
            if self.display_frame is None:
                self.display_frame = tk.Frame(self.root, bd=2, relief="groove")
                self.display_frame.grid(row=1, column=0, padx=10, pady=10)

                # Initialize image display with first molecule
                # Initialize sorting defaults
                self.sort_column = "CdId"
                self.sort_direction = "ASC"
                self.current_index = 0

                self.img_display = Label(self.display_frame)
                self.img_display.grid(row=1, column=0, padx=5, pady=10)

                self.info_display = tk.Frame(self.display_frame, bd=2, relief="groove")
                self.info_display.grid(row=1, column=1, padx=5, pady=10)

                self.refresh_display()

                
                self.database_path_label = Label(self.display_frame, text=db_path, border=2, relief="sunken")
                self.database_path_label.grid(row=0, column=1, padx=5, pady=10)

                self.next_img = Button(self.display_frame, text="Next Molecule", command=self.display_next)
                self.next_img.grid(row=2, column=2, padx=15, pady=10)

                icon_path = os.path.join(os.path.dirname(__file__), "..", "..", "images", "save-icon.png")
                img2 = Image.open(icon_path).resize((30, 30))
                self.save_icon = ImageTk.PhotoImage(img2)
                self.save_image = Button(self.display_frame, image=self.save_icon, command=self.save_current_image, width=30, height=30)
                self.save_image.grid(row=2, column=1, padx=15, pady=10)

                self.display_jump_label = Label(self.display_frame, text="Go to molecule #:", width=15)
                self.display_jump_label.grid(row=2, column=5, padx=5, pady=10)
                self.display_jump_entry = tk.Entry(self.display_frame, width=10, validate="key", validatecommand=(self.display_frame.register(self.validate_digit), "%P"))
                self.display_jump_entry.grid(row=2, column=4, padx=5, pady=10)
                self.display_jump_entry.bind("<Return>", lambda event: self.display_jump())

                self.chg_ord = Button(self.display_frame, text="Change Sort Order", command=self.update_order)
                self.chg_ord.grid(row=2, column=3, padx=15, pady=10)

                self.prev_img = Button(self.display_frame, text="Previous Molecule", command=self.display_previous)
                self.prev_img.grid(row=2, column=0, padx=15, pady=10)
                
                db_label = Label(self.display_frame, text="Database Path:")
                db_label.grid(row=0, column=0, padx=5, pady=10)

                options = ["CdId","Molecular weight", "LogP", "H-bond Donors", "H-bond Acceptors"]
                self.selected_option = tk.StringVar(self.display_frame)
                self.selected_option.set(options[0])
                dropdown = OptionMenu(self.display_frame, self.selected_option, *options, command=self.update_sort)
                dropdown.grid(row=0, column=2, padx=5, pady=10)
            else:
                # Update label if display frame already exists
                if self.database_path_label is not None:
                    self.database_path_label.config(text=db_path)


    def validate_digit(self, P):
        return P.isdigit() or P == ""
    
    def display_jump(self):
        if self.image_handler is None:
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            index = int(self.display_jump_entry.get()) - 1  # Convert to 0-based index
            if index < 0:
                messagebox.showerror("Error", "Please enter a positive integer!")
                return
            molecule_count = self.get_molecule_count()
            if index >= molecule_count:
                messagebox.showerror("Error", f"Please enter a number between 1 and {molecule_count}!")
                return
            self.current_index = index
            self.refresh_display()
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid integer!")

    def get_molecule_count(self):
        if self.database is None or self.database.db_file is None:
            messagebox.showerror("Error", "Database not initialized!")
            return 0
        
        con = sqlite3.connect(self.database.db_file)
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM molecules")
        count = cur.fetchone()[0]
        con.close()
        return count

if __name__ == "__main__":
    app()