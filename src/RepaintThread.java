public class RepaintThread extends Thread{
    Frame frame;
    RepaintThread(Frame frame){
        this.frame = frame;
    }

    @Override
    public void run() {
        while (true) {
            try {
                frame.repaint();
                if (frame.cameraXYZ[0] < 20){
                    frame.cameraXYZ[0] = frame.cameraXYZ[0] + 1;
                } else {
                    frame.cameraXYZ[0] = -20;
                }
                Thread.sleep(Frame.fpsMS); // Esperamos entre actualizaciones para obtener FPS
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
}
