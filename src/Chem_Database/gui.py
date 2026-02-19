import tkinter as tk
import os
from tkinter import filedialog, messagebox, Button, Label, OptionMenu
from Chem_Database.database import Database
from Chem_Database.image import ImageHandler
from Chem_Database.img_display import InfoHandler
from PIL import ImageTk, Image
import sys
import sqlite3
import logging

# Configure logging to use the same log file as database.py
logging.basicConfig(
    filename="database_errors.log",
    level=logging.WARNING,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def resource_path(relative_path):
    """
    Get the absolute path to a resource, works for development and for PyInstaller.
    
    _MEIPASS is an attribute added by PyInstaller to the sys module, which points to the temporary folder where the bundled application is extracted at runtime. This function checks if _MEIPASS exists and constructs the path accordingly, ensuring that resources are correctly located whether the code is run as a script or as a bundled executable.
    """
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)



class app:
    """
    Main application class for the Chem Database GUI.
    """
    def __init__(self):
        """
        Initialize the main application window and its components.
        """
        self.root = tk.Tk()
        self.root.title("Chem Database GUI")

        # Initialize instance variables
        self.file_path = None
        self.database = None
        self.image_handler = None
        self.display_frame = None
        self.database_path_label = None
        self.current_index = 1
        self.current_photo = None
        self.sort_column = "CdId" # Default sort column
        self.sort_direction = "ASC" # Default sort direction

        # Create file selection frame
        file_frame = tk.Frame(self.root, bd=2, relief="groove")
        file_frame.grid(row=0, column=0, padx=10, pady=10)

        '''SDF File selection components'''
        self.file_path_label = Label(file_frame, text="No file selected", border=2, relief="sunken")
        self.file_path_label.grid(row=0, column=1, padx=5, pady=10)

        file_select = Button(file_frame, text="Select SDF File", command=self.open_file)
        file_select.grid(row=0, column=0, padx=5, pady=10)

        load_sdf = Button(file_frame, text="Load SDF", command=self.load_sdf)
        load_sdf.grid(row=0, column=2, padx=5, pady=10)

        '''Database creation and loading components'''
        create_db = Button(file_frame, text="Create Database", command=self.create_database) # Button to create database from loaded SDF file and any applied filters
        create_db.grid(row=1, column=0, padx=5, pady=10)

        load_db = Button(file_frame, text="Load Database", command=self.load_database) # Button to load existing database, does not require SDF file to be loaded first and does not apply filters
        load_db.grid(row=1, column=2, padx=5, pady=10)
        

        '''Create collapsable details frame for filter options'''
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
        collapse_button = Button(file_frame, text="Filter Options", command=self.toggle_details) # Button to show/hide filter options, does not apply filters until "Create Database" button is clicked
        collapse_button.grid(row=1, column=1, padx=5, pady=10)

        self.root.mainloop() # Start the main event loop of the application
    
    def toggle_details(self):
        """
        Toggle the visibility of the filter options frame.
        """
        if self.details_frame.winfo_viewable():
            self.details_frame.grid_remove()
        else:
            self.details_frame.grid(row=3, column=0, columnspan=4, padx=10, pady=10)

    def open_file(self):
        """
        Open a file dialog to select an SDF file.
        """
        file_path = filedialog.askopenfilename(filetypes=[("SDF files", "*.sdf")])
        if file_path:
            self.file_path = file_path
            self.file_path_label.config(text=file_path)
            # Create a new Database instance with the selected file
            self.database = Database(self.file_path)
            messagebox.showinfo("File Selected", f"Selected: {file_path}")

    def load_sdf(self):
        """
        Load the selected SDF file into the database.
        """
        if self.file_path is None:
            logging.error("Load SDF failed: No SDF file selected")
            messagebox.showerror("Error", "Please select an SDF file first!")
            return
        
        if self.database is None:
            logging.error("Load SDF failed: Database not initialized")
            messagebox.showerror("Error", "Database not initialized!")
            return
        
        try:
            self.database.load_sdf()
            messagebox.showinfo("Success", "SDF file loaded successfully!")
        except Exception as e:
            logging.error(f"Failed to load SDF file: {str(e)}")
            messagebox.showerror("Error", f"Failed to load SDF file:\n{str(e)}")

    def create_database(self):
        """
        Create a database from the loaded SDF file with optional filters.
        """
        if self.file_path is None:
            logging.error("Create database failed: No SDF file selected")
            messagebox.showerror("Error", "Please select an SDF file first!")
            return
        
        if self.database is None or self.database.df is None:
            logging.error("Create database failed: SDF file not loaded")
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
            save_path = filedialog.asksaveasfilename(
                defaultextension=".db",
                filetypes=[("SQLite Database", "*.db")],
                initialfile="chem_database.db"
            ) # Prompt user to select save location for database, default name is "chem_database.db"

            if not save_path:
                return  # User cancelled

            self.database.db_file = save_path
            self.database.create_database(filter_list=filter_list) # Create database with applied filters, if any
            db_path = self.database.db_file # Get the path of the created database
            self.image_handler = ImageHandler(db_path) # Initialize ImageHandler with the created database to enable image retrieval for display
            # Create and show display frame after successful database creation
            if self.display_frame is None:
                self.initialize_database(db_path)
            else:
                # Update label if display frame already exists
                if self.database_path_label is not None:
                    self.database_path_label.config(text=db_path)
            
            messagebox.showinfo("Success", "Database created successfully!")
        except Exception as e:
            logging.error(f"Failed to create database: {str(e)}")
            messagebox.showerror("Error", f"Failed to create database:\n{str(e)}")
    
    def update_sort(self, selected_option):
        """
        Update the sort column based on the selected option.
        """
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
        """
        Toggle the sort direction between ascending and descending.
        """
        self.sort_direction = "DESC" if self.sort_direction == "ASC" else "ASC"
        self.current_index = 0  # Reset to first molecule when sort order changes
        self.refresh_display()
    
    def display_next(self):
        """
        Display the next molecule in the database based on the current sort order and index.
        """
        if self.image_handler is None:
            logging.error("Display next failed: Database not created")
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            # Increment to next molecule
            if self.current_index < self.get_molecule_count() - 1:
                self.current_index += 1
                self.refresh_display()
        except Exception as e:
            logging.error(f"Failed to display next image: {str(e)}")
            messagebox.showerror("Error", f"Failed to display image:\n{str(e)}")

    def display_previous(self):
        """
        Display the previous molecule in the database based on the current sort order and index.
        """
        if self.image_handler is None:
            logging.error("Display previous failed: Database not created")
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            # Decrement to previous molecule
            if self.current_index > 0:
                self.current_index -= 1
                self.refresh_display()
        except Exception as e:
            logging.error(f"Failed to display previous image: {str(e)}")
            messagebox.showerror("Error", f"Failed to display image:\n{str(e)}")
    
    def refresh_display(self):
        """
        Refresh the displayed molecule image and information based on the current index and sort order.
        """
        if self.image_handler is None:
            logging.error("Refresh display failed: Database not created")
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
        """
        Update the displayed information for the current molecule.
        """
        # clear previous info
        for widget in self.info_display.winfo_children():
            widget.destroy()
        
        if self.database is None or self.database.db_file is None:
            logging.error("Update info display failed: Database not initialized")
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
            logging.error(f"Failed to retrieve molecule information: {str(e)}")
            messagebox.showerror("Error", f"Failed to retrieve information:\n{str(e)}")
    
    def save_current_image(self):
        """
        Save the currently displayed molecule image to a file.
        """
        if self.image_handler is None:
            logging.error("Save image failed: Database not created")
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            img = self.image_handler.get_image_by_offset(self.current_index)
            save_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
            if save_path:
                img.save(save_path)
                messagebox.showinfo("Success", f"Image saved to {save_path}")
        except Exception as e:
            logging.error(f"Failed to save image: {str(e)}")
            messagebox.showerror("Error", f"Failed to save image:\n{str(e)}")
        
    def load_database(self):
        """
        Load a database from a file and initialize the image handler.
        """
        db_path = filedialog.askopenfilename(filetypes=[("SQLite Database", "*.db")])
        if db_path:
            self.database = Database(sdf_file=None, db_file=db_path)
            self.image_handler = ImageHandler(db_path)
            messagebox.showinfo("Database Loaded", f"Loaded database from: {db_path}")
            if self.display_frame is None:
                self.initialize_database(db_path)
            else:
                # Update label if display frame already exists
                if self.database_path_label is not None:
                    self.database_path_label.config(text=db_path)


    def validate_digit(self, P):
        """
        Validate that the input for the jump entry is a digit or empty.
        """
        return P.isdigit() or P == ""
    
    def display_jump(self):
        """
        Display the molecule at the specified index based on user input.
        """
        if self.image_handler is None:
            logging.error("Display jump failed: Database not created")
            messagebox.showerror("Error", "Database not created yet!")
            return
        
        try:
            index = int(self.display_jump_entry.get()) - 1  # Convert to 0-based index
            if index < 0:
                logging.error(f"Display jump failed: Invalid index {index + 1} (must be positive)")
                messagebox.showerror("Error", "Please enter a positive integer!")
                return
            molecule_count = self.get_molecule_count()
            if index >= molecule_count:
                logging.error(f"Display jump failed: Index {index + 1} out of range (max: {molecule_count})")
                messagebox.showerror("Error", f"Please enter a number between 1 and {molecule_count}!")
                return
            self.current_index = index
            self.refresh_display()
        except ValueError:
            logging.error(f"Display jump failed: Invalid integer input '{self.display_jump_entry.get()}'")
            messagebox.showerror("Error", "Please enter a valid integer!")

    def get_molecule_count(self):
        """
        Get the total number of molecules in the database.
        """
        if self.database is None or self.database.db_file is None:
            logging.error("Get molecule count failed: Database not initialized")
            messagebox.showerror("Error", "Database not initialized!")
            return 0
        
        con = sqlite3.connect(self.database.db_file)
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM molecules")
        count = cur.fetchone()[0]
        con.close()
        return count
    
    def initialize_database(self, db_path):
        """
        Initialize the database display frame and related widgets.
        """
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

        self.sorting_frame = tk.Frame(self.display_frame)
        self.sorting_frame.grid(row=0, column=0, columnspan=3, padx=5, pady=10)


        self.database_path_label = Label(self.sorting_frame, text=db_path, border=2, relief="sunken")
        self.database_path_label.grid(row=0, column=1, padx=5, pady=10)

        self.interact_frame = tk.Frame(self.display_frame)
        self.interact_frame.grid(row=2, column=0, columnspan=3, padx=5, pady=10)

        self.next_img = Button(self.interact_frame, text="Next Molecule", command=self.display_next) # Button to display the next molecule based on current sort order and index
        self.next_img.grid(row=0, column=4, padx=15, pady=10)

        save_icon_path = resource_path("images/save-icon.png")
        img2 = Image.open(save_icon_path).resize((30, 30))
        self.save_icon = ImageTk.PhotoImage(img2)
        self.save_image = Button(self.interact_frame, image=self.save_icon, command=self.save_current_image, width=30, height=30) # Button to save the currently displayed molecule image, uses a save icon image for the button
        self.save_image.grid(row=0, column=1, padx=15, pady=10)

        self.display_jump_label = Label(self.interact_frame, text="Go to molecule # :", width=15)
        self.display_jump_label.grid(row=0, column=2, padx=5, pady=10)
        self.display_jump_entry = tk.Entry(self.interact_frame, width=10, validate="key", validatecommand=(self.interact_frame.register(self.validate_digit), "%P")) # Entry widget for user to input molecule index to jump to, validates that input is a digit and triggers display_jump() on Enter key press
        self.display_jump_entry.grid(row=0, column=3, padx=5, pady=10)
        self.display_jump_entry.bind("<Return>", lambda event: self.display_jump()) # Bind Enter key to trigger display_jump() function when user inputs an index and presses Enter

        self.chg_ord = Button(self.sorting_frame, text="Change Sort Order", command=self.update_order) # Button to toggle the sort order between ascending and descending for the currently selected sort column, updates the displayed molecule accordingly
        self.chg_ord.grid(row=0, column=3, padx=15, pady=10)

        self.prev_img = Button(self.interact_frame, text="Previous Molecule", command=self.display_previous) # Button to display the previous molecule based on current sort order and index
        self.prev_img.grid(row=0, column=0, padx=15, pady=10)
        
        db_label = Label(self.sorting_frame, text="Database Path:")
        db_label.grid(row=0, column=0, padx=5, pady=10)

        options = ["CdId","Molecular weight", "LogP", "H-bond Donors", "H-bond Acceptors"]
        self.selected_option = tk.StringVar(self.sorting_frame)
        self.selected_option.set(options[0])
        dropdown = OptionMenu(self.sorting_frame, self.selected_option, *options, command=self.update_sort) # Dropdown menu to select the column to sort by, updates the displayed molecule based on the selected sort column and current sort direction
        dropdown.grid(row=0, column=2, padx=5, pady=10)

if __name__ == "__main__":
    app()