import os, sys, glob
import filecmp, shutil

def listFolder(fld):
    for f in glob.glob(fld+'/.*'):
        print f

def listFolderW(fld):
    for root, dirs, files in os.walk(fld):
        print root.replace(fld,'testcode'), len(files)

def compareFiles(orgF,newF):
    if os.path.exists(newF):
        if not filecmp.cmp(orgF,newF):
            print 'update', newF
    else:
        print 'new file', newF

def copyFiles(orgF,newF):
    if os.path.exists(newF):
        if not filecmp.cmp(orgF,newF):
            print 'update', newF
            shutil.copy(orgF,newF)
    else:
        print 'new file', newF
        shutil.copy(orgF,newF)
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
            if os.access(root+"/"+f, os.X_OK) and ('source' in root or 'validation' in root):
                continue
            if dryRun:
                compareFiles(root+"/"+f,newroot+'/'+f)
            else:
                copyFiles(root+"/"+f,newroot+'/'+f)
    if not dryRun:
        updateVersion(fld1)

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


def addVersion(fld1):
    currentDir=os.path.split(os.path.abspath(__file__))[0]
    os.chdir(fld1)
    os.system("git log -1 --format=commit' '%H > "+currentDir+"/bin/version")
    os.system("git describe --tags >> "+currentDir+"/bin/version")
    os.chdir(currentDir)
    os.system("svn add bin/version")


#listFolder('phosim/source')
#listFolderW('../../phosim/source')
#cpFolder('../../phosim/')
if len(sys.argv)<2:
    updateFolder('/home/epeng/phosim/')
else:
    if sys.argv[1] == 'False':
        updateFolder('/home/epeng/phosim/',dryRun=False)
    else:
        updateFolder('/home/epeng/phosim/')

#addVersion('../../phosim/')


#change source to src
# 1. configure
# 2. mv Makefile to src/.
# 3. update src/Makefile, validation/Makefile, validation/unittest.cpp
