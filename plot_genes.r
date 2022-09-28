
#gene plotting function from Raquel

plotGene<-function(start, end, strands, chroms,aLabel, aChrom,aMin, aMax, idCex = 0.8 , type = "gene", geneIdPatDif, labels = "top", useAxis = TRUE)
{
        myColours<-c("#0000ff","#ff0000")
        secondColours<-c("#ff0000")
        myOddColours<-
        aSelection<-which(chroms == aChrom & end > aMin & start < aMax)
        anXlim<-c(aMin, aMax)
        anYlim<-c(-0.5,5)
        if(type != "gene")
        {
                return("it only plots gene structures")
        }
        else
        {
                newStart<-start[aSelection]
                newEnd<-end[aSelection]
                newEnd<-newEnd[order(newStart)]
                newLabels<-aLabel[aSelection][order(newStart)]
                newStart<-sort(newStart)
                newStrand<-strands[aSelection]
                newSizes<-newEnd-newStart
#                idCex<-0.8
                suppressWarnings(par(mgp= c(1, 0.8, 0), mar= c(4, 1, 2, 1)))
                plot(1:10, type = "n",  ylab ="", xlab = "",frame = FALSE, main = "", xlim = anXlim,ylim =anYlim ,axes = FALSE)
		if (length(newStart) < 3)
                {
			yLocations<-rep(c(0.5, 0, 0, 0.5, 1.7, 1.2, 1.2, 1.7, 2.9,2.4,2.4,2.9, 4.1, 3.6, 3.6, 4.1))
                }
		if (length(newStart) >= 3)
                {
	                yLocations<-rep(c(0.5, 0, 0, 0.5, 1.7, 1.2, 1.2, 1.7, 2.9,2.4,2.4,2.9, 4.1, 3.6, 3.6, 4.1),round(length(newStart)/4))
		}
                toTakeHight<-0.25
                newColours<-NULL
                newColours[which(newStrand == "+" | newStrand == 1)]<-myColours[1]
                newColours[which(newStrand == "-" | newStrand == 0)]<-myColours[2]
                if(!is.null(geneIdPatDif))
                {
                        newColours[grep(geneIdPatDif, aLabel)]<-secondColours
                }
                midpoints<-newStart + newSizes/2
                range<-(aMax-aMin)/100
                labelsY<-NULL
                for ( i in 1:length(newStart))
                {
                        yLoc<-yLocations[((i-1)*4+1):(((i-1)*4)+4)]
                        if((newSizes[i])/2 <range )
                        {
                                newRange<-(newSizes[i])/2
                        }
                        else
                        {
                                newRange<-range
                        }
                        if(newStrand[i] == "-" | newStrand[i] == 0)
                        {
                                labelsY[i]<-yLoc[1]+0.2
                                polygon(c(newStart[i]+newRange, newStart[i], newStart[i]+newRange, rep(newEnd[i], each = 2)),c(yLoc[1],yLoc[1]-toTakeHight, yLoc[2:4]),  col = c(newColours[i]), border =newColours[i], lwd=1)
                        }
                        else
                        {
                                labelsY[i]<-yLoc[1]+0.2
                                polygon(c(rep(newStart[i], each = 2), newEnd[i]-newRange,newEnd[i],newEnd[i]-newRange ),c(yLoc[1:3],yLoc[4]-toTakeHight, yLoc[4]), col = c(newColours[i]), border =newColours[i], lwd=1)
                        }
                }
                if(labels == "top")
                {
                        text(midpoints, labelsY+0.2, newLabels, cex = idCex )
                }
                if(useAxis)
                {
                        megaNumber<-1000000
                        kbNumber<-1000
                        if(aMax > megaNumber)
                        {

                                axisStart<-round(aMin/megaNumber, 1 )
                                axisEnd<-round(aMax/megaNumber, 1)
                                axis(1,seq(axisStart*megaNumber,axisEnd*megaNumber,round((axisEnd*megaNumber-axisStart*megaNumber)/5, 1)),  lwd =1, lwd.ticks=1)
                        }
                        else if (aMax > kbNumber)
                        {
                                axisStart<-round(aMin/kbNumber, 1 )
                                axisEnd<-round(aMax/kbNumber, 1)
                                axis(1,seq(axisStart*kbNumber,axisEnd*kbNumber,round((axisEnd*kbNumber-axisStart*kbNumber)/5, 1)),  lwd =1, lwd.ticks=1)
                        }
                        else
                        {
                                axis(1,seq(aMin,aMax,round((aMax-aMin)/5, 1)), lwd =0, lwd.ticks=1)
                        }
                }
        }
}
