@echo On
@Rem delete all svn-related files
@PROMPT [Com]

@for /r . %%a in (.) do @if exist "%%a/.svn" rd /s /q "%%a/.svn"
@Rem for /r . %%a in (.) do @if exist "%%a/.svn" @echo "%%a/.svn"

@echo Mission Completed.
@pause