# MyProjects

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

A simple, secure client-server cloud storage system implemented in C using the **POSIX API**.

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

### Future Improvements:
- Replace fork() with threads for improved concurrency
- Fix issues related to large file transfers

---

## 3. Genetic Algorithms and Heuristic methods
The main problem approaced is finding the global minima of some computationally demanding functions, through experiments with genetic algorithms, simulated annealing algorithms and hill climbing algorithm.
    The results of the experiments are in the pdf file. The experiments for the Simulated Annealing method have been done with starting temperature=640 and sigma=400 and the Hill Climbing used 1000 repetitions.
