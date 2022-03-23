# STF, the simple trajectory format

This is part a rant, part a toy project.

## Design goals
In order of importance:

* **Simplicity** STF is designed with the main goal of being easy to read and write,
so readers and writers can be implemented in most programming languages.

* **Disk efficiency** It should be at least considerably better than DCD,
to make it worth it.

* **Write speed** Write operations should be reasonably fast.

* **Read speed**  Same for read operations.

### Why this order?

I think the single mayor factor in determining the "success" of a format, 
is the success of the program that uses it. A format might be a hot mess, but,
if the program that uses it is successful, so will the format. And the other way around.
Of course, nothing to do there. I actually implemented a variant of this format which
is compatible with the multi-XYZ format (not really _too_ different from the current
version) but the file size was pretty bad. Only a bit better than DCD, so not worthy.

Apart from that, I think simplicity is the main factor for a successful format. 
Of course this will not be a successful format, as in "other people will use it", but
simplicity allows me to quickly build readers for other languages (Python, particularly).
and make it compatible.

Disk space is also important. This is the problem with the DCD format.
