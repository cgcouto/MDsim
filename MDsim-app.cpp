#include "wx/sizer.h"
#include "wx/wx.h"

#include "MDsim-functions.h"

// this render loop is heavily based on code found here: https://wiki.wxwidgets.org/Making_a_render_loop

// store particle positions to be updated and rendered
long double ** currentParticles = importData("initial_particles.txt");

// scales the sim visualization to make it fit better on one screen
const double PICTURE_SCALE = 1.5;

class BasicDrawPane : public wxPanel
{
    
public:
    BasicDrawPane(wxFrame* parent);
	
    void paintEvent(wxPaintEvent& evt);
    void paintNow();
    void render( wxDC& dc );
    
    DECLARE_EVENT_TABLE()
};

class MyFrame;

class MyApp: public wxApp
{
    bool render_loop_on;
    bool OnInit();
    void onIdle(wxIdleEvent& evt);
    
    MyFrame* frame;
    BasicDrawPane* drawPane;
public:
    void activateRenderLoop(bool on);
        
};

IMPLEMENT_APP(MyApp)

class MyFrame : public wxFrame
{
public:
    MyFrame() : wxFrame((wxFrame *)NULL, -1,  wxT("MDSimApp"), wxPoint(50,50), 
                        wxSize(WIDTH/PICTURE_SCALE,HEIGHT/PICTURE_SCALE))
    {
    }
    void onClose(wxCloseEvent& evt)
    {
        wxGetApp().activateRenderLoop(false);
        evt.Skip(); // don't stop event, we still want window to close
    }
    DECLARE_EVENT_TABLE()
};


BEGIN_EVENT_TABLE(MyFrame, wxFrame)
EVT_CLOSE(MyFrame::onClose)
END_EVENT_TABLE()

bool MyApp::OnInit()
{
    render_loop_on = false;
    
    wxBoxSizer* sizer = new wxBoxSizer(wxHORIZONTAL);
    frame = new MyFrame();
	
    drawPane = new BasicDrawPane( frame );
    sizer->Add(drawPane, 1, wxEXPAND);
	
    frame->SetSizer(sizer);
    frame->Show();
    
    activateRenderLoop(true);
    return true;
} 

void MyApp::activateRenderLoop(bool on)
{
    if(on && !render_loop_on)
    {
        Connect( wxID_ANY, wxEVT_IDLE, wxIdleEventHandler(MyApp::onIdle) );
        render_loop_on = true;
    }
    else if(!on && render_loop_on)
    {
        Disconnect( wxEVT_IDLE, wxIdleEventHandler(MyApp::onIdle) );
        render_loop_on = false;
    }
}
void MyApp::onIdle(wxIdleEvent& evt)
{
    if(render_loop_on)
    {
        drawPane->paintNow();
        evt.RequestMore(); // render continuously, not only once on idle
    }
}


BEGIN_EVENT_TABLE(BasicDrawPane, wxPanel)
EVT_PAINT(BasicDrawPane::paintEvent)
END_EVENT_TABLE()



BasicDrawPane::BasicDrawPane(wxFrame* parent) :
wxPanel(parent)
{
}


void BasicDrawPane::paintEvent(wxPaintEvent& evt)
{
    wxPaintDC dc(this);
    render(dc);
}

void BasicDrawPane::paintNow()
{
    wxClientDC dc(this);
    render(dc);
}

void BasicDrawPane::render( wxDC& dc )
{
    currentParticles = doOneFrame(currentParticles); // update the particle positions
    dc.SetBrush(*wxBLUE_BRUSH);

    // render each particle 
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        dc.DrawCircle( wxPoint(currentParticles[i][0]/PICTURE_SCALE,currentParticles[i][1]/PICTURE_SCALE), 
                      (PARTICLE_DIAM/2)/PICTURE_SCALE);
    }
}