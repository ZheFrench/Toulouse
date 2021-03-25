## How to use git


'''shell 

# Create git repository, the remote and push
git init  .
git remote add origin your_repo.git 
# The files in your dir are tracked only. 
git add .
# Files are committed, means any changes will be followed now.
git commit -m "First commit"
git push origin

#List the files and their states.

git status -s

#Add only modified files.

git add -u

# Give you the list of remote repository
git remote -v

# Create you upstream dir
git remote add origin https://github.com/ZheFrench/Toulouse.git

# Rename the url if you need...
git remote set-url origin https://github.com/ZheFrench/Toulouse.git

# Send code
git push origin
#id:ZheFrench/pwd:Camille34

# If you previously delete the file not using git rm , the file is staying in the remote...do that and i will disapear during your next push
git rm --cached 'file name'

# You did a mistake and push
git revert --no-commit <commit>

'''
