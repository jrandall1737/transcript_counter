# Transcript Counter --- Getting Started Guide

This project uses a Jupyter Notebook to analyze transcript files and
generate summary statistics.

Follow the steps below to set up your environment and run the notebook.

---

## Prerequisites

Before you begin, make sure you have:

- Python 3.8 or newer installed
- Jupyter Notebook or JupyterLab installed
- `pip` available in your environment

(Optional but recommended: use a virtual environment.)

---

## 1. Clone or Download the Project

If using Git:

```bash
git clone https://github.com/jrandall1737/transcript_counter
cd transcript_counter
```

Or download and extract the ZIP file into a local directory.

---

## 2. Set Up a Virtual Environment (Recommended)

Create and activate a virtual environment:

### Windows

```bash
python -m venv venv
venv\Scripts\activate
```

### macOS / Linux

```bash
python3 -m venv venv
source venv/bin/activate
```

---

## 3. Install Dependencies

Install all required Python packages from `requirements.txt`:

```bash
pip install -r requirements.txt
```

Verify installation:

```bash
pip list
```

---

## 4. Configure File Paths

Open the notebook in VsCode with the jupyter plugin installed:

Then open:

    transcript_counter.ipynb

Inside the notebook, locate the code section labeled:

    Setup file paths

Edit any file or directory paths in this section so they match the
locations on your system.

Make sure all referenced paths exist and are accessible.

---

## 5. Run the Notebook

Once dependencies are installed and file paths are configured:

1.  In VsCode, open `transcript_counter.ipynb`
2.  Select **Run All**,\
    or run each cell manually from top to bottom

The notebook will process the transcripts and generate the output files.

## Troubleshooting

### Package Installation Errors

Try upgrading pip:

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

---

### Jupyter Not Found

Install Jupyter if needed:

```bash
pip install notebook
```

or

```bash
pip install jupyterlab
```

---

### File Not Found Errors

- Double-check paths in **Setup file paths**
- Use absolute paths if relative paths cause issues
- Confirm permissions on target folders

---

## Support

If you encounter issues:

- Verify your Python version
- Confirm all dependencies installed correctly
- Recheck file paths
- Review error messages in notebook output

For additional help, contact the project maintainer.
