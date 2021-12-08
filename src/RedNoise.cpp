#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>

#define WIDTH 320
#define HEIGHT 240

using namespace std;
using namespace glm;

vec3 camera(0.0,0.0,4.0)

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}
//nov numberofvalues
std::vector<CanvasPoint> interpolate(CanvasPoint from, CanvasPoint to, int nov){
	std::vector<CanvasPoint> cp;
	float xNov = (to.x - from.x)/(nov - 1);
	float yNov = (to.y - from.y)/(nov - 1);
	CanvasPoint v;

	for(float i = 0; i < nov -1; i++){
		v.x = v.x + x_nov;
		v.y = v.y + y_nov;
		cp.push_back(v);
	}
	return cp;


}

void drawLine(CanvasPoint from, CanvasPoint to, glm::vec3 colour, DrawingWindow& window ) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float nos =  std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/nos;
	float yStepSize = yDiff/nos;

	uint32_t col = (255 << 24) + ((int)colour[0] << 16) + ((int)colour[1] << 8) + (int)colour[2];
	for (float i=0.0; i<nos; i++){
		float x = from.x + (xStepSize * i);
		float y = from.y + (yStepSize * i);
		window.setPixelColour(round(x), round(y), col)

	}

}



void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
