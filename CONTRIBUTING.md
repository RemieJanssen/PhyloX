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

Then make some changes and push the branch.
```
git add [your changed file]
git commit -m "[some message describing the changes]"
git push -u origin [yourname]-[description_of_branch]
```

Now create a merge request for your branch. Make sure to set `RemieJanssen/PhyloX:main` as the base, i.e. the branch to merge into. See [this page](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork) on the GitHub docs for more guidance.
