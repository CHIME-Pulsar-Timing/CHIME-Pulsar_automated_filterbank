The python presto functions are really useful, but a lot have lines like "from presto import x" in them.
I want to be able to run them without having to install presto (& there are many machines where I can't)
AND a lot of the python presto modules actually only depend on each other rather than the core presto routines

So I'm gathering them here and tweaking them to avoid presto imports
(or replace them with presto_without_presto)
