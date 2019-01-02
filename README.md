# RenewableSimulator

The renewable energy simulator

## Known Issues

#### ssc.dylib library not found on MacOS

[Reference](http://log.zyxar.com/blog/2012/03/10/install-name-on-os-x/)

After successful compilation, if you receive the following error message when you try to run the executable, you probably need to change the install name of your `ssc.dylib` file.

```
$ ./RenewableSimulator 
dyld: Library not loaded: ssc.dylib
  Referenced from: /Users/wuh20/github/RenewableSimulator/build/./RenewableSimulator
  Reason: image not found
Abort trap: 6
```


Let's check whether the install name is not correctly set up for the `ssc.dylib`. In my case, I put this library under `/Users/wuh20/packages/sam/osx64`.

```
$ otool -D ssc.dylib 
ssc.dylib:
ssc.dylib
```

In this case, the install name of the library is just itself, not the full path to the file. This is causing problems during the runtime because the program that is calling this library cannot find this library. Let's change its install name to the full path.

```
$ install_name_tool -id /Users/wuh20/packages/sam/osx64/ssc.dylib ssc.dylib 
$ otool -D ssc.dylib 
ssc.dylib:
/Users/wuh20/packages/sam/osx64/ssc.dylib
$ otool -L ssc.dylib 
ssc.dylib:
	/Users/wuh20/packages/sam/osx64/ssc.dylib (compatibility version 0.0.0, current version 0.0.0)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1252.0.0)
	/usr/lib/libc++.1.dylib (compatibility version 1.0.0, current version 400.9.0)
```

As you can see, after we change the install name of the library, the first line of the output of shared library used is changed accordingly. Now you should go back and rerun your program.
