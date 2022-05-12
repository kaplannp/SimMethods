import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation

class DataLoader:

    def load(self, fileName):
        '''
        parse a file and return three 3d np arrays, i, j, p.
        shape (n_runs, nRows, nCols). (arr[0] is a slice of one run)
        '''
        with open(fileName, 'r') as file:
            shape = self.getArrShape(file)
            iArr = np.empty(shape)
            jArr = np.empty(shape)
            pArr = np.empty(shape)
            file.seek(0)
            file.readline() #don't want first sep
            fileString = file.read();
            blocks = fileString.split('---\n')
            for i, block in enumerate(blocks):
                if(i%3==0):
                    iArr[i//3] = self.parseBlock(block)
                elif(i%3==1):
                    jArr[i//3] = self.parseBlock(block)
                else: #(i%3==2):
                    pArr[(i//3)] = self.parseBlock(block)
            return (iArr, jArr, pArr)


    def parseBlock(self, block):
        '''
        takes as input a string representing a single grid single run.
        loads that string into numpy
        @param block single grid run
        @returns array representing the grid shape i, j
        '''
        lines = block[:-1].split('\n') #drop the last newline character
        shapedList = [line.split('  ') for line in lines]
        return np.array(shapedList)

    def countRuns(self, file):
        '''
        @param file an open file object for reading
        @return count for the number of runs in the file
        Postcondition: file seeker has been modified
        '''
        file.seek(0);
        fileString = file.read();
        count = fileString.count("---\n")
        assert (count % 3 == 0),"invalid number of --- in file"
        return int(count/3)

    def getArrShape(self, file):
        '''
        @param file an open file object for reading
        @return a shape for arrays
        Postcondition: file seeker has been modified
        '''
        nRuns = self.countRuns(file)
        file.seek(0);
        file.readline();
        line = file.readline();
        nCols = len(line.split('  '))
        nRows = 0;
        while(line != "---\n"):
            nRows+=1;
            line = file.readline();
        return (nRuns, nRows, nCols)


class Plotter:

    def plotFieldMovie(self, iField, jField):
        fig = plt.figure()
        anim = FuncAnimation(fig, self.plotField,
                range(iField.shape[0]), fargs=(iField, jField))
        return anim

    def plotField(self, i, iField, jField):
        plt.clf();
        plt.quiver(jField[i], iField[i], scale = 1e2)




loader = DataLoader();
iArr, jArr, pArr = loader.load("simData.out")
plotter = Plotter()
while True:
    anim = plotter.plotFieldMovie(iArr, jArr)
    #writergif = animation.PillowWriter(fps=5) 
    #anim.save("withReflection.gif", writergif)
    plt.show()

