# Chemical Structure Database Creator and GUI
### Current Version:
V1.0.0

## About
This program creates a GUI to take an SDF file with optional filters, which is used to populate a database; each entry is then displayed individually with calculated stats. 

## How To Use
There are multiple ways to run this program, the easiest being to download the latest release and run the executable file. 
You can also either:
- Download the repo from github and run either gui.py or main.py directly
- Clone the repo and install via pip, then run by simply typing 'chem-database-tool' in the terminal while in the correct working directory. Exact instructions can be seen blow

### Installing via PIP
1. Clone the repo:
```bash
git clone https://github.com/BeyondOurMinds/Chem-Comp-Database.git
```

2. Switch to correct directory:
```bash
cd Chem-Comp-Database
```

3. Create a virtual environment:
```bash
python -m venv .venv
```

4. Activate virtual environment:

**Windows (PowerShell):**
```powershell
.venv\Scripts\Activate.ps1
```

**Windows (Command Prompt):**
```cmd
.venv\Scripts\activate.bat
```

**Linux/Mac:**
```bash
source .venv/bin/activate
```

5. Install the package:
```bash
pip install .
```

6. Run the package:
```bash
chem-database-tool
```