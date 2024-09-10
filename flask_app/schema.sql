DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id INTEGER PRIMARY KEY AUTOINCREMENT,
    project_name TEXT UNIQUE NOT NULL,
    project_author TEXT NOT NULL,
    creation_date DATETIME NOT NULL,
    directory_path TEXT NOT NULL,
    vdj_path TEXT NOT NULL,
    adata_path TEXT NOT NULL
);