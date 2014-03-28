# Content

Readme for VIC source code formating using [uncrustify](http://uncrustify.sourceforge.net/), which can also be installed through most package managers (e.g. [Macports](http://www.macports.org)).

# Purpose 

To create easy to read, constistently formatted source code and avoid unnecessary commits because of formatting and whitespace changes.

# Usage
To format a VIC source code file:

    uncrustify -c uncrustify_VIC_c.cfg [files ...]

There are other command line options available:

    man uncrustify

To format all code in a directory without creating a backup (assuming your shell is `bash`):

    for file in *[ch]
    do
        echo $file
        uncrustify -c uncrustify_VIC_c.cfg --no-backup --replace $file
    done
    
# Guidelines

To make life easy:

* **do** run uncrustify on the file[s] that are part of a pull request (before you commit)

* **don't** change the configuration file

* **don't** change the formatting after you run uncrustify

If you really don't like the VIC settings, you can create your own configuration file for uncrustify and run the files through that before you edit them. Just remember to run them through the VIC version before you commit (and don't commit whitespace changes).

