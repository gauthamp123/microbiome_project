import sqlite3

conn = sqlite3.connect("hmmtop.db")

c = conn.cursor()
c.execute("SELECT name FROM sqlite_master WHERE type='table'")
print(c.fetchall())