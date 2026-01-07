`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 01/07/2026 09:57:12 AM
// Design Name: 
// Module Name: Demo_tb_Hough_Transform
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module Demo_tb_Hough_Transform ();

    parameter IMG_WIDTH     = 640;
    parameter IMG_HEIGHT    = 480;
    parameter IO_DATA_WIDTH = 24;
    parameter FRAME_SIZE    = IMG_WIDTH * IMG_HEIGHT;
    

    reg aresetn;
    reg ap_clk;
    

    reg [IO_DATA_WIDTH-1:0] tx_data_reg;
    reg                     tx_valid_reg;
    reg                     tx_last_reg;
    reg                     tx_user_reg;

    wire [IO_DATA_WIDTH-1:0] s_axis_tdata;  assign s_axis_tdata  = tx_data_reg;
    wire                     s_axis_tvalid; assign s_axis_tvalid = tx_valid_reg;
    wire                     s_axis_tlast;  assign s_axis_tlast  = tx_last_reg;
    wire                     s_axis_tuser;  assign s_axis_tuser  = tx_user_reg;
    
    wire                    s_axis_tready; 


    wire [IO_DATA_WIDTH-1:0] m_axis_tdata;
    wire                     m_axis_tvalid; 
    reg                      m_axis_tready; 
    wire                     m_axis_tlast;
    wire                     m_axis_tuser;
    

    integer i, j;
    reg [IO_DATA_WIDTH-1:0] input_image [0:FRAME_SIZE-1];
    integer output_file;

    
    

    Demo_Hough_Transform_Design_Source #(
        .IMG_WIDTH(IMG_WIDTH),
        .IMG_HEIGHT(IMG_HEIGHT)
        
    ) dut (
        .aresetn(aresetn),
        .ap_clk(ap_clk),
        .s_axis_tdata(s_axis_tdata),
        .s_axis_tvalid(s_axis_tvalid),
        .s_axis_tready(s_axis_tready),
        .s_axis_tlast(s_axis_tlast),
        .s_axis_tuser(s_axis_tuser),
        .m_axis_tdata(m_axis_tdata),
        .m_axis_tvalid(m_axis_tvalid),
        .m_axis_tready(m_axis_tready),
        .m_axis_tlast(m_axis_tlast),
        .m_axis_tuser(m_axis_tuser)
    );

localparam INPUT_FILENAME = "input_image_Hough.hex";
localparam OUTPUT_FILENAME = "output_image_Hough.hex";

    initial begin
        ap_clk = 0;
        forever #5 ap_clk = ~ap_clk; 
    end


    initial begin
      
        $readmemh(INPUT_FILENAME, input_image);
        output_file = $fopen(OUTPUT_FILENAME, "w");

        
        aresetn       = 0;     
        tx_valid_reg  = 0;
        tx_last_reg   = 0;
        tx_user_reg   = 0;
        tx_data_reg   = 0;
        m_axis_tready = 1;    

      
        #2; 
        @(negedge ap_clk); 

        repeat (10) @(posedge ap_clk); 
        
        @(posedge ap_clk); 
        aresetn = 1;       

        $display("TESTBENCH: Starting simulation. Reset deasserted at %t ns.", $time);
        
        fork
           
            begin
               
                tx_valid_reg = 1;
                tx_user_reg  = 1; 
                tx_data_reg  = 0; 
                tx_last_reg  = 0; 

                @(posedge ap_clk); 
                

                tx_user_reg  = 0; 
                tx_valid_reg = 0; 

                $display("TESTBENCH: Sent SOF trigger at %t ns. Waiting for DUT to finish reset...", $time);
                

                wait (s_axis_tready === 1'b1); 
                @(posedge ap_clk); 
                $display("TESTBENCH: DUT is ready for frame input at %t ns. Starting pixel data transmission.", $time);
                

                for (i = 0; i < FRAME_SIZE; i = i + 1) begin
                  
                    wait (s_axis_tready === 1'b1);
                    
                    tx_valid_reg = 1;
                    tx_data_reg  = input_image[i];

                    if ((i + 1) % IMG_WIDTH == 0) begin
                        tx_last_reg = 1;
                    end else begin
                        tx_last_reg = 0;
                    end
                    
                    @(posedge ap_clk); 
                end
                
                tx_valid_reg = 0;
                tx_last_reg  = 0;

                $display("TESTBENCH: Finished sending input frame at time %t ns.", $time);
            end

            begin

                for (j = 0; j < FRAME_SIZE; j = j + 1) begin
                    wait (m_axis_tvalid === 1'b1); 
                    
                    if (j == 0) begin
                        $display("TESTBENCH: Received first output pixel at time %t ns.", $time);
                    end

                    $fdisplay(output_file, "%h", m_axis_tdata);
                    
                    @(posedge ap_clk); 
                end
                
                m_axis_tready = 0; 
                $display("TESTBENCH: Successfully received all output pixels at %t ns.", $time);
            end
        join 
        
        $fclose(output_file);
        $display("TESTBENCH: Test finished successfully at time %t ns.", $time);
        #100; 
        $finish;
    end
    
endmodule
