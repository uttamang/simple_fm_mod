For Windoze
How to use git and github
- Download and install https://git-for-windows.github.io/
- Open cmd
- go to a directory (empty is best i guess) with cd
- do "git init"
- do "git remote add origin https://github.com/stoertebeker23/simple_fm_mod"
- do "git pull origin" -> git dir should be in the selected folder now

When you change a file
- do "git status"  and check if the file is marked as modified
- do "git add <file>"
- do "git commit -m "<message about the work>" "
- do "git push origin HEAD"
- keep commits small

Proxy
- do "git config --global http.proxy http://<username>:<pw>@proxy.hs-karlsruhe.de:8888"
- do the same for https
