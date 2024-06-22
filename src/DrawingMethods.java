import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.List;
import java.util.Queue;
import java.util.*;

public class DrawingMethods {
    int width, height;
    public Color originalColor;
    Frame frame;
    BufferedImage buffer;

    public DrawingMethods(Frame frame, BufferedImage buffer) {
        this.frame = frame;
        this.buffer = buffer;
        this.width = buffer.getWidth();
        this.height = buffer.getHeight();
        originalColor = Color.getColor(String.valueOf(buffer.getRGB(0, 0)));
    }

    public void putPixel(int x, int y, Color a) {
        if (x >= 0 && x < width && y >= 0 && y < height) {
            // Establece el color del píxel en el buffer
            buffer.setRGB(x, y, a.getRGB());
        }
    }

    //......................Lines......................

    public void horizontalLine(int x1, int x2, int y, Color a){
        int dx = x2 - x1;
        int leftToRight = (int) Math.signum(dx);
        dx = Math.abs(dx);
        int x = x1;
        for (int i = 0; i < dx; i++){
            x += leftToRight;
            putPixel(x, y, a);
        }
    }

    public void verticalLine(int y1, int y2, int x, Color a){
        int dy = y2 - y1;
        int topToBottom = (int) Math.signum(dy);
        dy = Math.abs(dy);
        int y = y1;
        for (int i = 0; i < dy; i++){
            y += topToBottom;
            putPixel(x, y, a);
        }
    }

    public void bresenham(int x1, int y1, int x2, int y2, Color a) {
        int dx = (x2 - x1);
        int dy = (y2 - y1);
        int leftToRight = (int) Math.signum(dx);
        int topToBottom = (int) Math.signum(dy);
        int doubleDx = 2 * Math.abs(dx);
        int doubleDy = 2 * Math.abs(dy);
        float m = (float) dy / dx;
        int x = x1;
        int y = y1;
        int p;
        if (Math.abs(m) > 1){
            p = doubleDx - Math.abs(dy);
            for (int i = 0; i < Math.abs(dy); i++) {
                putPixel(x, y, a);
                if (p <= 0) {
                    p += doubleDx;
                } else {
                    p += doubleDx - doubleDy;
                    x += leftToRight;
                }
                y += topToBottom;
            }
        } else {
            p = doubleDy - Math.abs(dx);
            for (int i = 0; i < Math.abs(dx); i++) {
                putPixel(x, y, a);
                if (p <= 0) {
                    p += doubleDy;
                } else {
                    p += doubleDy - doubleDx;
                    y += topToBottom;
                }
                x += leftToRight;
            }
        }
    }

    public void line(int[][] xy, Color color){
        for (int i = 0; i < xy[0].length - 1; i++){
            bresenham(xy[0][i],xy[1][i],xy[0][i+1],xy[1][i+1],color);
        }
    }

    //......................Shapes......................

    public void rectangle(int x1, int y1, int x2, int y2, Color color, boolean fill){
        verticalLine(y1,y2,x1,color);
        horizontalLine(x1,x2,y2,color);
        verticalLine(y2,y1,x2,color);
        horizontalLine(x2,x1,y1,color);
        if (fill){
            double[][] polygonPoints = {
                    {x1, y1},
                    {x2, y1},
                    {x2, y2},
                    {x1, y2}
            };
            scanLine(polygonPoints,color);
        }
    }

    public void polygon(int numPoints, int[] xCoordinates, int[] yCoordinates, Color a, boolean fill) {
        if (xCoordinates.length != yCoordinates.length) {
            System.out.println("Los arreglos de coordenadas deben tener la misma cardinalidad");
            return;
        }
        int x1, x2, y1, y2;
        if (fill) {
            for (int i = 0; i < numPoints; i++) {
                x1 = xCoordinates[i];
                y1 = yCoordinates[i];

                if ((i + 1) < numPoints) {
                    x2 = xCoordinates[(i + 1)];
                    y2 = yCoordinates[(i + 1)];
                } else {
                    x2 = xCoordinates[0];
                    y2 = yCoordinates[0];
                }
                bresenham(x1, x2, y1, y2, a);
            }
            int rows = xCoordinates.length;

            double[][] polygonPoints = new double[rows][2]; // 2 columnas

            for (int i = 0; i < rows; i++) {
                polygonPoints[i][0] = xCoordinates[i]; // Llenar la primera columna
                polygonPoints[i][1] = yCoordinates[i]; // Llenar la segunda columna
            }

            scanLine(polygonPoints,a);
        }
        else {
            for (int i = 0; i < numPoints; i++) {
                x1 = xCoordinates[i];
                y1 = yCoordinates[i];

                if ((i + 1) < numPoints) {
                    x2 = xCoordinates[(i + 1)];
                    y2 = yCoordinates[(i + 1)];
                } else {
                    x2 = xCoordinates[0];
                    y2 = yCoordinates[0];
                }

                //midpoint(x1, x2, y1, y2, a);
                bresenham(x1, x2, y1, y2, a);
                //dda(x1, x2, y1, y2, a);
            }
        }
    }

    public void circleMidPoint(int xc, int yc, int r, Color a, boolean fill) {
        int x = 0;
        int y = r;
        int d = 1 - r;
        int deltaE = 3;
        int deltaSE = -2 * r + 5;
        putPixel(xc + x, yc + y, a);
        putPixel(xc - x, yc + y, a);
        putPixel(xc + x, yc - y, a);
        putPixel(xc - x, yc - y, a);
        putPixel(xc + y, yc + x, a);
        putPixel(xc - y, yc + x, a);
        putPixel(xc + y, yc - x, a);
        putPixel(xc - y, yc - x, a);

        while (y > x) {
            if (d < 0) {
                d += deltaE;
                deltaE += 2;
                deltaSE += 2;
            } else {
                d += deltaSE;
                deltaE += 2;
                deltaSE += 4;
                y--;
            }
            x++;
            drawCirclePoints(xc, yc, x, y, a);
        }

        if (fill){
            //flood(xc,yc,a);
            floodFillIterative(xc,yc,a,r,r,xc,yc);
        }
    }

    private void drawCirclePoints(int xc, int yc, int x, int y, Color a) {
        putPixel(xc + x, yc + y, a);
        putPixel(xc - x, yc + y, a);
        putPixel(xc + x, yc - y, a);
        putPixel(xc - x, yc - y, a);
        putPixel(xc + y, yc + x, a);
        putPixel(xc - y, yc + x, a);
        putPixel(xc + y, yc - x, a);
        putPixel(xc - y, yc - x, a);
    }

    public void ellipse(int xc, int yc, int rx, int ry, double angle, Color color, boolean fill, boolean precision){
        int estimatedPerimeter = (int) (2 * Math.PI * Math.max(rx, ry));
        while (angle >= 360) angle -= 360;
        double alpha = Math.toRadians(angle);
        int x, y;
        for (int i = 0; i < estimatedPerimeter; i++) {
            double j = (i * 2 * Math.PI) / estimatedPerimeter;
            x = xc + (int) ((double) rx * Math.cos(j) * Math.cos(alpha) - (double) ry * Math.sin(j) * Math.sin(alpha));
            y = yc + (int) ((double) rx * Math.cos(j) * Math.sin(alpha) + (double) ry * Math.sin(j) * Math.cos(alpha));
            putPixel(x,y,color);
        }
        if (fill){
            if (rx > 40 || ry > 40) {
                scanLineEllipse(xc,yc,rx,ry,color);
            } else {
                if (!precision){
                    flood(xc,yc,color);
                } else {
                    floodFillIterative(xc, yc, color, rx, ry, xc, yc);
                }
            }
        }
    }

    //......................Fill Algorithms......................

    public void scanLine(double[][] polygonPoints, Color fillColor) {
        double minY = Double.MAX_VALUE;
        double maxY = Double.MIN_VALUE;
        for (double[] point : polygonPoints) {
            minY = Math.min(minY, point[1]);
            maxY = Math.max(maxY, point[1]);
        }

        for (int y = (int) minY; y <= maxY; y++) {
            List<Double> xIntersections = new ArrayList<>();

            for (int i = 0; i < polygonPoints.length; i++) {
                double[] startPoint = polygonPoints[i];
                double[] endPoint = polygonPoints[(i + 1) % polygonPoints.length];

                if ((startPoint[1] <= y && endPoint[1] >= y) || (endPoint[1] <= y && startPoint[1] >= y)) {
                    double dx = endPoint[0] - startPoint[0];
                    double dy = endPoint[1] - startPoint[1];
                    if (dy != 0) {
                        double x = startPoint[0] + (y - startPoint[1]) * dx / dy;
                        xIntersections.add(x);
                    }
                }
            }

            Collections.sort(xIntersections);

            for (int i = 0; i < xIntersections.size() - 1; i += 2) {
                int startX = xIntersections.get(i).intValue();
                int endX = xIntersections.get(i + 1).intValue();
                horizontalLine(startX, endX, y, fillColor);
            }
        }
    }

    public void scanLineEllipse(int xc, int yc, int rx, int ry, Color fillColor) {
        for (int y = yc - ry; y <= yc + ry; y++) {
            for (int x = xc - rx; x <= xc + rx; x++) {
                if (isPointInEllipse(x, y, rx, ry, xc, yc)) {
                    putPixel(x, y, fillColor);
                }
            }
        }
    }

    public void flood(int x, int y, Color newColor) {
        int width = buffer.getWidth();
        int height = buffer.getHeight();

        Queue<Integer> queue = new ArrayDeque<>();
        queue.add(x);
        queue.add(y);

        while (!queue.isEmpty()) {
            int curX = queue.poll();
            int curY = queue.poll();

            if (curX < 0 || curX >= width || curY < 0 || curY >= height) {
                continue;
            }

            int pixelColor = buffer.getRGB(curX, curY);
            if (pixelColor == newColor.getRGB() ){
                continue;
            }

            putPixel(curX, curY, newColor);

            queue.add(curX + 1);
            queue.add(curY);
            queue.add(curX - 1);
            queue.add(curY);
            queue.add(curX);
            queue.add(curY + 1);
            queue.add(curX);
            queue.add(curY - 1);
        }
    }

    public void floodSameColor(int x, int y, Color newColor) {
        int originalColor = buffer.getRGB(x,y);
        int width = buffer.getWidth();
        int height = buffer.getHeight();

        Queue<Integer> queue = new ArrayDeque<>();
        queue.add(x);
        queue.add(y);

        while (!queue.isEmpty()) {
            int curX = queue.poll();
            int curY = queue.poll();

            if (curX < 0 || curX >= width || curY < 0 || curY >= height) {
                continue;
            }

            int pixelColor = buffer.getRGB(curX, curY);
            if (pixelColor != originalColor || pixelColor == newColor.getRGB()){
                continue;
            }

            putPixel(curX, curY, newColor);

            queue.add(curX + 1);
            queue.add(curY);
            queue.add(curX - 1);
            queue.add(curY);
            queue.add(curX);
            queue.add(curY + 1);
            queue.add(curX);
            queue.add(curY - 1);
        }
    }

    public void floodFillIterative(int x, int y, Color fillColor, int a, int b, int centerX, int centerY) {
        Set<Point> visitedPixels = new HashSet<>();

        if (!isPointInEllipse(x, y, a, b, centerX, centerY)) {
            return;
        }

        Queue<Point> queue = new LinkedList<>();
        queue.add(new Point(x, y));

        while (!queue.isEmpty()) {
            Point p = queue.remove();
            x = p.x;
            y = p.y;

            // Verificar si el píxel ya ha sido visitado
            if (visitedPixels.contains(p) || !isPointInEllipse(x, y, a, b, centerX, centerY)) {
                continue;
            }

            // Marcar el píxel como visitado
            visitedPixels.add(p);

            // Dibujar el píxel
            putPixel(x, y, fillColor);

            // Agregar píxeles adyacentes a la cola si aún no han sido visitados
            queue.add(new Point(x + 1, y));
            queue.add(new Point(x - 1, y));
            queue.add(new Point(x, y + 1));
            queue.add(new Point(x, y - 1));
        }
    }

    private boolean isPointInEllipse(int x, int y, int a, int b, int centerX, int centerY) {
        double dx = x - centerX;
        double dy = y - centerY;
        return (dx * dx) / (a * a) + (dy * dy) / (b * b) <= 1;
    }

    //......................Curves......................

    public void sinSinglePeriod(int xi, int yi, double amplitude, double wavelenght, int steps, double scale, Color color){
        int[][] points = new int[steps+1][2];
        int c = 0;
        double y;
        double step = Math.PI / ((double) steps);
        for (double i = 0; Math.abs(i) <= Math.PI; i += step){
            y = yi - (Math.sin(i * 1/wavelenght) / amplitude) * scale;
            points[c][0] = (int) Math.round(i * scale) + xi;
            points[c][1] = (int) Math.round(y);
            c++;
        }
        line(points,color);
    }

    public int[][] twoPointArc(int x1, int y1, int x2, int y2, int ry, int points) {
        int xc = (int) (Math.min(x1, x2) + Math.abs(x2 - x1) / (double) 2);
        int yc = (int) (Math.min(y1, y2) + Math.abs(y2 - y1) / (double) 2);

        int rx = (int) (Math.sqrt(Math.pow(Math.abs(x2 - x1), 2) + Math.pow(Math.abs(y2 - y1), 2))) / 2;
        if (points == 0) {
            points = (int) (Math.PI * Math.max(rx, ry));
        }

        int[][] matrix = new int[2][points];
        double dx = x2 - x1;
        double dy = y2 - y1;
        double alpha = dx != 0 ? Math.atan(dy / dx) : Math.PI / 2;

        if (x1 > x2 || (x1 == x2 && y1 > y2)) ry *= -1;

        for (int i = 0; i < points; i++) {
            double t = (i * Math.PI) / points;

            int x = xc + (int) ((double) rx * Math.cos(t) * Math.cos(alpha) - (double) ry * Math.sin(t) * Math.sin(alpha));
            int y = yc + (int) ((double) rx * Math.cos(t) * Math.sin(alpha) + (double) ry * Math.sin(t) * Math.cos(alpha));

            matrix[0][i] = x;
            matrix[1][i] = y;
        }
        return matrix;
    }

    //......................Transformations......................//

    //--------Rotations--------//

    public void polygonRotation(int xPivot, int yPivot, int[][] coordinates, double angle, Color color, boolean fill){
        int points = coordinates[0].length;
        while (angle >= 360) angle -= 360;
        int[][] matrixResult = new int[2][points];

        angle = 0 + Math.toRadians(angle);

        for (int i = 0; i < points; i++){
            matrixResult[0][i] = xPivot + (int) (Math.round(coordinates[0][i]-xPivot)*Math.cos(angle)-(coordinates[1][i]-yPivot)*Math.sin(angle));
            matrixResult[1][i] = yPivot + (int) (Math.round(coordinates[0][i]-xPivot)*Math.sin(angle)+(coordinates[1][i]-yPivot)*Math.cos(angle));
        }

        line(matrixResult,color);
    }

    public void pointRotation(int xPivot, int yPivot, int x, int y, double angle, Color color, boolean fill){
        int newX,newY;
        while (angle >= 360) angle -= 360;
        angle = 0 + Math.toRadians(angle);

        newX = xPivot + (int) (Math.round(x-xPivot)*Math.cos(angle)-(y-yPivot)*Math.sin(angle));
        newY = yPivot + (int) (Math.round(x-xPivot)*Math.sin(angle)+(y-yPivot)*Math.cos(angle));

        if (fill){
            flood(newX,newY,color);
        }
    }

    public int[][] getRotatedPoint(int xPivot, int yPivot, int x, int y, double angle){
        int newX,newY;
        int[][] point = new int[2][1];

        while (angle >= 360) angle -= 360;
        angle = 0 + Math.toRadians(angle);

        newX = xPivot + (int) (Math.round(x-xPivot)*Math.cos(angle)-(y-yPivot)*Math.sin(angle));
        newY = yPivot + (int) (Math.round(x-xPivot)*Math.sin(angle)+(y-yPivot)*Math.cos(angle));

        point[0][0] = newX;
        point[1][0] = newY;
        return point;
    }

    public void ellipseRotation(int xPivot, int yPivot, int xc, int yc, int rx, int ry, double angle, Color color, boolean fill, boolean precitionFill){
        int newX,newY;
        while (angle >= 360) angle -= 360;
        double angleRad = 0 + Math.toRadians(angle);

        newX = xPivot + (int) (Math.round(xc-xPivot)*Math.cos(angleRad)-(yc-yPivot)*Math.sin(angleRad));
        newY = yPivot + (int) (Math.round(xc-xPivot)*Math.sin(angleRad)+(yc-yPivot)*Math.cos(angleRad));

        ellipse(newX,newY,rx,ry,angle,color,fill,precitionFill);
    }

    //--------Scalations--------//

    public void polygonScalation(int xPivot, int yPivot, int[][] coordinates, double sx, double sy, Color color, boolean fill){
        int points = coordinates[0].length;
        int[][] matrixResult = new int[2][points];

        for (int i = 0; i < points; i++){
            matrixResult[0][i] = xPivot + (int) (Math.round((coordinates[0][i] - xPivot)*sx));
            matrixResult[1][i] = yPivot + (int) (Math.round((coordinates[1][i] - yPivot)*sy));
        }

        line(matrixResult,color);
    }

    public void pointScalation(int xPivot, int yPivot, int x, int y, double sx, double sy, Color color, boolean fill){
        int newX,newY;

        newX = xPivot + (int) (Math.round((x - xPivot)*sx));
        newY = yPivot + (int) (Math.round((y - yPivot)*sy));

        if (fill){
            flood(newX,newY,color);
        }
    }

    public int[][] getScaledPoint(int xPivot, int yPivot, int x, int y, double sx, double sy){
        int newX,newY;
        int[][] point = new int[2][1];

        newX = xPivot + (int) (Math.round((x - xPivot)*sx));
        newY = yPivot + (int) (Math.round((y - yPivot)*sy));

        point[0][0] = newX;
        point[1][0] = newY;
        return point;
    }

    public void ellipseScalation(int xPivot, int yPivot, int xc, int yc, int rx, int ry, double sx, double sy, double angle, Color color, boolean fill,boolean precision){
        int newX,newY;
        int newRx = (int) Math.round(rx * sx);
        int newRy = (int) Math.round(ry * sy);

        newX = xPivot + (int) (Math.round((xc - xPivot)*sx));
        newY = yPivot + (int) (Math.round((yc - yPivot)*sy));

        ellipse(newX,newY,newRx,newRy,angle,color,fill,precision);

    }

    //--------Escalation & rotation--------//

    public void polygonScalationRotation(int xPivot, int yPivot, int[][] coordinates, double sx, double sy, double angle, Color color, boolean fill){
        int points = coordinates[0].length;
        int[][] matrixResult = new int[2][points];

        for (int i = 0; i < points; i++){
            matrixResult[0][i] = xPivot + (int) (Math.round((coordinates[0][i] - xPivot)*sx));
            matrixResult[1][i] = yPivot + (int) (Math.round((coordinates[1][i] - yPivot)*sy));
        }

        polygonRotation(xPivot,yPivot,matrixResult,angle,color,fill);
    }

    public void pointScalationRotation(int xPivot, int yPivot, int x, int y, double sx, double sy, Color color, boolean fill){
        int newX,newY;

        newX = xPivot + (int) (Math.round((x - xPivot)*sx));
        newY = yPivot + (int) (Math.round((y - yPivot)*sy));

        if (fill){
            flood(newX,newY,color);
        }
    }

    public void ellipseScalationRotation(int xPivot, int yPivot, int xc, int yc, int rx, int ry, double sx, double sy, double angle, Color color, boolean fill,boolean precitionFill){
        int newX,newY;
        int newRx = (int) Math.round(rx * sx);
        int newRy = (int) Math.round(ry * sy);

        newX = xPivot + (int) (Math.round((xc - xPivot)*sx));
        newY = yPivot + (int) (Math.round((yc - yPivot)*sy));

        ellipse(newX,newY,newRx,newRy,angle,color,fill,precitionFill);
    }
    public void prismProjection(int[][] prism, double xp, double yp, double zp, int xPivot, int yPivot, double scale, Color color){
        double[] u = new double[prism[0].length];
        int vertices = (prism[0].length)/2;
        System.out.println(vertices);
        int[][] projectedPrism = new int[2][prism[0].length];
        for (int i = 0; i < prism[0].length; i++){
            u[i] = (double) (-prism[2][i]) /zp;
            projectedPrism[0][i] = xPivot + (int) Math.round((prism[0][i] + xp*u[i])* scale);
            projectedPrism[1][i] = yPivot - (int) Math.round((prism[1][i] + yp*u[i])* scale);
        }
        for (int i = 0; i < vertices; i++){
            if (i < vertices - 1){
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][i+1],projectedPrism[1][i+1],color);
                bresenham(projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],projectedPrism[0][i+vertices+1],projectedPrism[1][i+vertices+1],color);
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],color);
            } else {
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][0],projectedPrism[1][0],color);
                bresenham(projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],projectedPrism[0][vertices],projectedPrism[1][vertices],color);
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],color);
            }
        }
    }

    public void prismPerspective(int[][] prism, double xp, double yp, double zp, int xPivot, int yPivot, double scale, Color color){
        double[] u = new double[prism[0].length];
        int vertices = (prism[0].length)/2;
        System.out.println(vertices);
        int[][] projectedPrism = new int[2][prism[0].length];
        for (int i = 0; i < prism[0].length; i++){
            u[i] = -zp / (prism[2][i] - zp);
            projectedPrism[0][i] = xPivot + (int) Math.round((xp + (prism[0][i] - xp)*u[i]) * scale);
            projectedPrism[1][i] = yPivot - (int) Math.round((yp + (prism[1][i] - yp)*u[i]) * scale);
        }
        for (int i = 0; i < vertices; i++){
            if (i < vertices - 1){
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][i+1],projectedPrism[1][i+1],color);
                bresenham(projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],projectedPrism[0][i+vertices+1],projectedPrism[1][i+vertices+1],color);
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],color);
            } else {
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][0],projectedPrism[1][0],color);
                bresenham(projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],projectedPrism[0][vertices],projectedPrism[1][vertices],color);
                bresenham(projectedPrism[0][i],projectedPrism[1][i],projectedPrism[0][i+vertices],projectedPrism[1][i+vertices],color);
            }
        }
    }

}
