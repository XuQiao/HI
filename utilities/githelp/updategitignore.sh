#!/bin/sh
#To untrack every file that is now in your .gitignore:

#First commit any outstanding code changes, and then, run this command:

git rm -r --cached .
#This removes any changed files from the index(staging area), then just run:

git add .
#Commit it:

git commit -m ".gitignore is now working"
