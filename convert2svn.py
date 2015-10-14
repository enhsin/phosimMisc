import os, sys, glob
import filecmp, shutil
import difflib

def listFolder(fld):
    for f in glob.glob(fld+'/.*'):
        print f

def listFolderW(fld):
    for root, dirs, files in os.walk(fld):
        print root.replace(fld,'testcode'), len(files)

def compareFiles(oldF,newF):
    if os.path.exists(newF):
        if not filecmp.cmp(oldF,newF):
            print 'update', newF
    else:
        print 'new file', newF, oldF

def copyFiles(oldF,newF):
    if os.path.exists(newF):
        if not filecmp.cmp(oldF,newF):
            print 'update', newF
            shutil.copy(oldF,newF)
    else:
        print 'new file', newF
        shutil.copy(oldF,newF)
        os.system("svn add "+newF)

def updateFolder(fld1,dryRun=True):
    for root, dirs, files in os.walk(fld1):
        if 'git' in root or 'condor' in root:
            continue

        if 'source' in root:
            newroot=root.replace(fld1,'').replace('source','src/source')
        else:
            newroot=root.replace(fld1,'')

        if newroot == '':
            for f in files:
                if f in ['configure','Makefile','Doxyfile']:
                    continue
                elif '.csh' in f or '.lis' in f:
                    continue
                elif f in ['COPYING','README']:
                    newroot=''
                else:
                    newroot='src/tobin/'
                if dryRun:
                    compareFiles(root+"/"+f,newroot+f)
                else:
                    copyFiles(root+"/"+f,newroot+f)

            continue
        elif newroot == 'validation':
            for f in files:
                if f in ['.gitignore','Makefile','unittest.o','unittest']:
                    continue
                elif f in ['unittest.h','unittest.cpp']:
                    newroot='src/validation'
                else:
                    newroot='validation'
                if dryRun:
                    compareFiles(root+"/"+f,newroot+'/'+f)
                else:
                    copyFiles(root+"/"+f,newroot+'/'+f)

            continue
        if 'docs' in root: #don't copy 2.2G docs
            continue
        if 'work' in root or 'output' in root or 'tool' in root:
            continue
        if 'bin' in root:
            continue

        if not os.path.exists(newroot):
            print 'no', newroot
            if not dryRun:
                os.system("mkdir "+newroot)
                os.system("svn add "+newroot)

        for f in files:
            if 'git' in f:
                continue
            if 'akefile' in f:
                continue
            if '.o' in f or '.pyc' in f:
                continue
            if '.tar' in f:
                continue
            if '.csh' in f:
                continue
            if os.access(root+"/"+f, os.X_OK) and ('source' in root or 'validation' in root) and 'cpp' not in f:
                continue
            if dryRun:
                compareFiles(root+"/"+f,newroot+'/'+f)
            else:
                copyFiles(root+"/"+f,newroot+'/'+f)
    if not dryRun:
        updateMakefile(currentPhosim)
        updateVersion(fld1)
    else:
        compareMakefile(currentPhosim)

def cpFolder(fld1):
    for root, dirs, files in os.walk(fld1):
        if 'git' in root:
            print root
            continue
        newroot=root.replace(fld1,'')
        if root == '../../phosim/':
            for f in files:
                os.system("cp "+root+"/"+f+" "+f)
                os.system("svn add "+f)
            continue
        if root == '../../phosim/docs': #don't copy 2.2G docs
            continue

        if not os.path.exists(newroot):
            print newroot
            os.system("mkdir "+newroot)
            os.system("svn add "+newroot)
        for f in files:
            if 'git' in f:
                continue
            os.system("cp "+root+"/"+f+" "+newroot+"/"+f)
            os.system("svn add "+newroot+"/"+f)

    addVersion(fld1)


def updateVersion(fld1):
    currentDir=os.path.split(os.path.abspath(__file__))[0]
    os.chdir(fld1)
    os.system("git log -1 --format=commit' '%H > "+currentDir+"/src/tobin/version")
    os.system("git describe --tags >> "+currentDir+"/src/tobin/version")

def addVersion(fld1):
    currentDir=os.path.split(os.path.abspath(__file__))[0]
    os.chdir(fld1)
    os.system("git log -1 --format=commit' '%H > "+currentDir+"/bin/version")
    os.system("git describe --tags >> "+currentDir+"/bin/version")
    os.chdir(currentDir)
    os.system("svn add bin/version")

def compareMakefile(fld1):
    dirs = [f for f in os.listdir(fld1+'source/') if '.tar' not in f and '.git' not in f]
    for dir in dirs:
        oldMF=fld1+'source/'+dir+'/makefile'
        if not os.path.exists(oldMF):
            oldMF=fld1+'source/'+dir+'/Makefile'
        newMF='src/source/'+dir+'/makefile'
        if not os.path.exists(newMF):
            newMF='src/source/'+dir+'/Makefile'
        oldF=open(oldMF).readlines()
        newF=open(newMF).readlines()
        result=difflib.unified_diff(oldF,newF)
        resultStr=''.join(list(result))
        print resultStr

def updateMakefile(fld1):
    dirs = [f for f in os.listdir(fld1+'source/') if '.tar' not in f and '.git' not in f]
    for dir in dirs:
        oldMF=fld1+'source/'+dir+'/makefile'
        if not os.path.exists(oldMF):
            oldMF=fld1+'source/'+dir+'/Makefile'
        newMF='src/source/'+dir+'/makefile'
        if not os.path.exists(newMF):
            newMF='src/source/'+dir+'/Makefile'
        f=open(newMF,"r+")
        newd=f.readlines()
        f.seek(0)
        i=0
        for line in open(oldMF).readlines():
            if 'cp ' in line:
                continue
            if 'CFLAGS' in line and 'CC' not in line:
                while 'CFLAGS' not in newd[i]:
                    i+=1
                f.write(newd[i])
                continue
            if 'LFLAGS' in line and 'CC' not in line:
                while 'LFLAGS' not in newd[i]:
                    i+=1
                f.write(newd[i])
                continue
            f.write(line)
        f.truncate()
        f.close()


currentPhosim=sys.argv[1]
if not os.path.exists(currentPhosim):
    print currentPhosim, 'does not exist'
    sys.exit()
if currentPhosim[-1]!='/':
    currentPhosim+='/'

if len(sys.argv)<3:
    updateFolder(currentPhosim)
else:
    if sys.argv[2] == 'False':
        updateFolder(currentPhosim,dryRun=False)
    else:
        updateFolder(currentPhosim)

