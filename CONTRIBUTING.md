We welcome any contributions, especially from researchers in mathematical phylogenetics.
As the project has just started and currently consists mostly of the scripts from one PhD project, we have no strict guidelines yet.

## Rough guideline
If you have ideas for improving the python package, simply open a new issue.
You may then also submit a pull request with the proposed changes.
Please link these to eachother by mentioning one in the other.

### Issues
An issue should consist of a proposed change, a reason for that change, and preferably a solution.
A solution is by no means required, but without it, it can take a lot longer for the change to take effect.

### Pull request
First fork the repo by clicking the "Fork" button on https://github.com/RemieJanssen/PhyloX. Then clone the main branch of your fork and checkout a new branch.
```
git clone git@github.com:<your_github_username>/PhyloX.git
cd phylox
git checkout -b [yourname]-[description_of_branch]
```
Then make some changes and push the branch. In the resulting output, there should be a link to github for opening a pull request for your newly pushed branch. Make sure to set `RemieJanssen:main` as the base, i.e. the branch to merge into.
```
git add [your changed file]
git commit -m "[some message describing the changes]"
git push -u origin [yourname]-[description_of_branch]
```
