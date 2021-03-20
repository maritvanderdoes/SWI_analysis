# How to use GitHub

## Starting with GitHub
The first time, you should clone the repository from the website with the following command: <code> git clone (link of the repository) </code>

## Commands
Following times you can:
- Check which branch you are, and the status of the branch: <code>git status </code>
- Download the recent github version of the code: <code> git fetch </code> or <code> git pull </code>
- Check which branches are there:  <code> git branch --list </code>
- create branch and go there <code>git branch (name_branch) -> git checkout (name_branch)</code>  or  <code>git checkout -b (name_branch)</code>
- When you are done with working on the branch, save the file, go back to the terminal and check the status of the folder <code>git status </code>
- Adjusted files are in red. Add the changes:  <code> git add (name_of_the_file) </code>
- When new feature is implemented, commit the changes with a message what is improved: <code> git commit -m "commit message" </code>
- Upload the newest version of your branch to gitub : <code> git push origin (namebranch) </code>
- When you are finished with the issue, merge your branch to the main branch on github by doing a pull request
- Delete your branch  <code> git branch -d (branchname) </code>


## Other useful links
tutorial:
https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet

list with commands
https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html

https://gist.github.com/rxaviers/7360908
