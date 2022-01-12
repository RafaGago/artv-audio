Basically from this StackOverflow answer:

https://stackoverflow.com/questions/39938253/how-to-properly-debug-a-cross-compiled-windows-code-on-linux

Install the cross debugger client and server:
> apt install gdb-mingw-w64 gdb-mingw-w64-target

Launch the server through Wine:
> wine Z:/usr/share/win64/gdbserver.exe localhost:1234 ~/.wine/drive_c/Program\ Files/REAPER\ \(x64\)/reaper.exe

Connect through GDB:
> x86_64-w64-mingw32-gdb C:>/Program\ Files/REAPER\ \(x64\)/reaper.exe


JSON for VSCode:
```
        {
            "type": "gdb",
            "request": "attach",
            "name": "Wine-reaper",
            "executable": "home/s0001192/.wine/drive_c/Program\ Files/REAPER\ \(x64\)/reaper.exe ",
            "target": ":1234",
            "remote": true,
            "cwd": "${workspaceRoot}",
            "valuesFormatting": "parseText",
            "gdbpath": "/usr/bin/x86_64-w64-mingw32-gdb",
            "setupCommands": [{
                "description": "path for standard libraries",
                "text": "-gdb-set solib-search-path ${workspaceRoot}/<TODO>",
                "ignoreFailures": true
            }]

        }
```

Notice that building with xwin on Debug mode is not possible, as Microsoft doesn't redistribute e.g "msvcrtd.dll".
