from numpy import *
from sys import argv

def calculateGraphDetails(linkMatrix):
    # Number of pages considered 
    numPages = len(linkMatrix)

    # Matrix to hold incoming link matrix for each node
    iLinks = [[] for i in range(numPages)]

    # Number of outgoing links in each node
    nLinksPerPage = zeros(numPages, int32)

    # Array to contain dangling nodes
    leafNodes = []

    # Calculate the sparse matrix
    for i in range(numPages):
        if(len(linkMatrix[i]) == 0):
            leafNodes.append(i)
        else:
            nLinksPerPage[i] = len(linkMatrix[i])
            for j in linkMatrix[i]:
                iLinks[j].append(i)

    iLinks = [array(i) for i in iLinks]
    nLinksPerPage = array(nLinksPerPage)
    leafNodes = array(leafNodes)

    return iLinks, nLinksPerPage, leafNodes
    
# Page rank calculation 
def pageRankCalc(incomingLinkArray, pageLinksArray, leafNodesArray, alpha, convergence, steps):
    nPages = len(incomingLinkArray)
    lNodes = leafNodesArray.shape[0]

    prNew = ones((nPages,), float32) / nPages
    prOld = ones((nPages,), float32) / nPages

    flag = False

    while not flag:
        for s in range(steps):
            prOld, prNew = prNew, prOld

            # Eigen value 1 x I vector 
            eigenValue = (1 - alpha) * sum(prOld) / nPages

            # Hyper link matrix
            hyperLinkValue = 0.0
            if lNodes > 0:
                hyperLinkValue = alpha * sum(prOld.take(leafNodesArray, axis = 0)) / nPages


            i = 0
            while i < nPages:
                page = incomingLinkArray[i]
                h = 0
                if page.shape[0]:
                    h = alpha * dot(
                        prOld.take(page, axis = 0),
                        1. / pageLinksArray.take(page, axis = 0)
                        )
                # PR value for each of the nodes
                prNew[i] = h + eigenValue + hyperLinkValue
                i = i + 1

        difference = prNew - prOld
        flag = (sqrt(dot(difference, difference)) / nPages < convergence)

        yield prNew


def pageRank(linkMatrix):
    alpha = 0.85
    convergence = 0.01
    steps = 10
    # Generate incoming link array as the transpose of the given input array of outgoing link matrix
    # Calculate the number of links in each of the nodes. 
    # leafNodesArray contains the nodes that do not have outgoing link
    incomingLinkArray, pageLinksArray, leafNodesArray = calculateGraphDetails(linkMatrix)


    # calculate the page rank for the input
    for pr in pageRankCalc(incomingLinkArray, pageLinksArray, leafNodesArray, alpha, convergence, steps):
        finalValue = pr
    return finalValue

def main():
    # User input
    filename = argv[1]
    inputText = open(filename, 'r')
    inputStr = inputText.read()

    inputStr = inputStr.split("\n")
    linkMatrix = []

    for line in inputStr:
        if (inputStr.index(line) == 0):
            numNodes = line.split(" ")[0]
            numEdges = line.split(" ")[1]
            for n in range(int(numNodes)):
                linkMatrix.append([])
        elif(line != ''):
            i = int(line.split(" ")[0])
            j = int(line.split(" ")[1])
            linkMatrix[i-1].append(j-1)
        else:
            print "The result is"
    print pageRank(linkMatrix)



# Start the execution
main()