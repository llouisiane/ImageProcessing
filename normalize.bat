::windows only: drop arbitrary files on the icon of this .bat file to create histogram normalised pngs
for %%x in (%*) do (
    "bin\Release\normalize.exe" %%x
)
pause