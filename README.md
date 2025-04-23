# MyProjects

Welcome to **MyProjects**, a collection of small-to-medium sized applications that explore topics such as 3D modeling, cloud storage, and heuristic optimization algorithms. Each project is designed to showcase practical implementations in areas like systems programming, computer graphics, and artificial intelligence.

---

## 1. 3D Viewer

A lightweight application for viewing and manipulating simple 3D objects in real time.

### 🎮 Controls:
- `↑` (Up Arrow) — Zoom in  
- `↓` (Down Arrow) — Zoom out  
- `Right Click + Drag` — Rotate camera (change point of view)  
- `Left Click + Drag`  
  - On object center: Move entire object  
  - On vertex: Move selected vertex  
- `Middle Click + Drag` — Rotate object around its origin (when clicked at center)

---

## 2. Cloud Storage

A simple, secure **client-server cloud storage** system implemented in C using the **POSIX API**.

### Features:
- C (Linux)
- Socket programming
- File handling
- POSIX system calls
- User login and session management
- File upload/download functionalities

### Supported Client Commands:
```bash
login               # Log into an existing account
cont nou            # Create a new account
exit                # Exit the client
logout              # Log out of the current session
upload              # Upload a file to the cloud
download <filename> # Download a file from the cloud
list                # View all stored files
delfile <filename>  # Delete a file
```

---

## 3. Genetic Algorithms and Heuristic methods
The main problem approaced is finding the global minima of some computationally demanding functions, through experiments with genetic algorithms, simulated annealing algorithms and hill climbing algorithm.
    The results of the experiments are in the pdf file.
