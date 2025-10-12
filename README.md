# 🌀 Christoffel Symbol Calculator

A symbolic calculator built with **SymPy** and **Streamlit** to compute Christoffel symbols for arbitrary spacetime metrics in General Relativity.

## Features

- Compute Christoffel symbols from any metric tensor
- Predefined GR metrics: Schwarzschild, Kerr, FLRW, Minkowski, etc.
- Supports custom metric input
- Symbolic simplification and LaTeX output
- Easy web UI via Streamlit

## Check out the **Christoffel Symbol Calculator** here:
[Christoffel Symbol Calculator](https://christoffel-symbol-calculator-fzqbkxctohbhcnvvnjnpfz.streamlit.app/)

## Project structure 
christoffel-symbol-calculator/
- app.py                - Streamlit UI for user input and output display
- christoffel.py        - Core logic to compute Christoffel symbols
- requirements.txt      - Python dependencies
- README.md             - This file

## To run this app locally :
- Clone this repository
git clone https://github.com/karmakar-rahul
/christoffel-symbol-calculator.git

- Navigate to the project directory
cd christoffel-symbol-calculator

- (Optional) Create a virtual environment
python -m venv venv
venv\Scripts\activate     # On Windows
source venv/bin/activate  # On Linux/Mac

- Install dependencies
pip install -r requirements.txt

- Run the app
streamlit run app.py

## Technologies used : 
- Python 3.x
- Streamlit
- SymPy

## Author
Rahul Karmakar, MSc Physics (Astrophysics Specialization)

## If you find this tool helpful, give it a star on GitHub! It helps others discover it too :)



