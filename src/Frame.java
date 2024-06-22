import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;

public class Frame extends JFrame {
    private BufferedImage mainBuffer;
    private DrawingMethods dmMain;
    private BufferedImage backgroundBuffer;
    private DrawingMethods dmBackground;
    int width = 800;
    int height = 800;
    int gridSpacing = 20;
    int[] cameraXYZ = {5,5,5};
    public static int fps;
    public static int fpsMS;

    Graphics2D gMain;

    private RepaintThread repaintThread;

    public Frame() {
        setTitle("Pixel");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        fps = 24;
        fpsMS = Math.round((float) 1000 /fps);

        setSize(width, height);

        setLayout(null);
        setLocationRelativeTo(null);

        setVisible(true);
        repaintThread = new RepaintThread(this);
        repaintThread.start();
    }

    @Override
    public void paint(Graphics g) {
        if (backgroundBuffer == null) {
            backgroundBuffer = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_ARGB);
            dmBackground = new DrawingMethods(this, backgroundBuffer);
            Graphics2D gBackground = backgroundBuffer.createGraphics();
            gBackground.setClip(0, 0, getWidth(), getHeight());
            gBackground.setColor(Color.gray);
            gBackground.fillRect(0,0,getWidth(),getHeight());
            drawCartesianPlane(gBackground, width, height, gridSpacing);
            int[][] cube = {
                    {0,0,5,5,0,0,5,5},
                    {0,5,5,0,0,5,5,0},
                    {0,0,0,0,5,5,5,5}
            };
            int[][] tetrisPrice = {
                    {0,0,5,5,10,10,5,5,-5,-5,0,0,5,5,0,0},
                    {0,30,30,20,20,10,10,0,0,30,30,20,20,10,10,0},
                    {0,0,5,5,10,10,5,5,5,5,10,10,15,15,10,10}
            };
            int[][] hexagonalPyramid = {
                    {-10,0,2,2,0,0,2,2},
                    {0,2,2,0,0,2,2,0},
                    {0,0,0,0,3,3,3,3}
            };
            gBackground.dispose();
        }
        if (mainBuffer == null){
            mainBuffer = new BufferedImage(getWidth(),getHeight(),BufferedImage.TYPE_INT_ARGB);
            dmMain = new DrawingMethods(this, mainBuffer);
            gMain = mainBuffer.createGraphics();
        }
        updateBuffer(g);
    }

    private void updateBuffer(Graphics g) {
        int[][] cubePerspective = {
                {0,0,2,2,0,0,2,2},
                {0,2,2,0,0,2,2,0},
                {0,0,0,0,3,3,3,3}
        };
        gMain.drawImage(backgroundBuffer,0,0,this);
        dmMain.prismPerspective(cubePerspective,cameraXYZ[0],cameraXYZ[1],cameraXYZ[2],400,400,20,Color.ORANGE);
        g.drawImage(mainBuffer, 0, 0, this);
    }

    public static void drawCartesianPlane(Graphics g, int width, int height, int gridSpacing) {
        g.setColor(Color.LIGHT_GRAY); // Líneas de la cuadrícula en gris claro

        // Dibujar líneas verticales
        for (int x = 0; x < width; x += gridSpacing) {
            g.drawLine(x, 0, x, height);
        }

        // Dibujar líneas horizontales
        for (int y = 0; y < height; y += gridSpacing) {
            g.drawLine(0, y, width, y);
        }

        g.setColor(Color.black); // Ejes en negro

        // Dibujar eje X
        g.drawLine(0, height / 2, width, height / 2);

        // Dibujar eje Y
        g.drawLine(width / 2, 0, width / 2, height);
    }
}
