# tumor-growth-graphs
# Tumor Growth Modeling - Logistic Equation Analysis

This project simulates tumor development using the **Logistic Growth Model**, a fundamental equation in mathematical biology for representing population growth under resource constraints. The project includes both analytical and numerical (RK45) solutions, parameter sensitivity analysis, and a chemotherapy intervention simulation.

---

## 📊 Model Overview

The growth of a tumor is modeled using the logistic differential equation, which accounts for the carrying capacity of the environment:

$$\frac{dP}{dt} = r \cdot P \cdot \left(1 - \frac{P}{K}\right)$$

**Key Parameters:**
* **$P(t)$**: Tumor size at time $t$ ($mm^3$).
* **$r$**: Proliferation rate (growth speed).
* **$K$**: Carrying capacity (maximum size).

---

## 🚀 Key Features

The simulation consists of four distinct analytical perspectives:
* **Analytical vs. Numerical Solution**: Validation of the `SciPy` `solve_ivp` (RK45) method against the exact analytical solution.
* **Effect of Growth Rates ($r$)**: Comparison of how different growth speeds reach carrying capacity.
* **Initial Size ($P_0$) Impact**: Analysis of how the starting volume affects growth timing.
* **Chemotherapy Simulation**: A dynamic simulation where $r$ is reduced on Day 10 to model treatment.

---

## 📈 Visualizations

![Tumor Growth Analysis](tumor_buyumesi_grafikler.jpg)

* **Top Left**: Displays the high accuracy between analytical and numerical models.
* **Top Right**: Shows that higher $r$ values lead to faster saturation.
* **Bottom Left**: Demonstrates that larger initial tumors hit the rapid growth phase earlier.
* **Bottom Right**: Visualizes the "treatment effect" showing a flattening of the curve after Day 10.

---

## 💻 Setup and Usage

### Prerequisites
* Python 3.x
* NumPy, Matplotlib, SciPy

### Running the Simulation
```bash
python tumor.py
